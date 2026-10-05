"""Check the off-screen lymphoid organs against the kinetics they are built from.

1. Thymic involution reproduces the measured TREC decline (>95% from age 25 to 60).
2. Anatomical distance to the draining node is NOT rate-limiting - transit is dominated
   by interstitial crawling, not path length.
3. The priming loop runs on measured timings: DC transit, 3-phase priming, expansion.
4. Cognate precursor frequency is indication-appropriate: HeLa is allogeneic (~7%),
   the real indications are autologous neoantigens (~1e-5, Alanio et al. 2010).
5. Marrow granulopoiesis matches the measured human production rate.
"""
import numpy as np

from tme.organs import (BODY_MASS, MARROW, NODE_BASIN, NODE_DISTANCE_CM, PRIMING,
                        THYMUS, thymic_export_per_day)
from tme.tcells import EFFECTOR
from tme.world import World

print("1. Thymic involution (measured: TRECs fall >95% between 25 and 60)")
y25, y60 = thymic_export_per_day(25), thymic_export_per_day(60)
for a in (20, 25, 40, 60, 80):
    print(f"   age {a:3d}: {thymic_export_per_day(a):.2e} cells/day")
print(f"   decline 25 -> 60: {100 * (1 - y60 / y25):.1f}%  (target >95%)")
print(f"   NOTE: at 60 the naive pool is maintained by peripheral division, not this output "
      f"(den Braber 2012)")

print("\n2. Is distance to the draining node rate-limiting?")
for ind in ("PDAC", "LUAD", "OV", "BRCA"):
    w = World(seed=1, indication=ind)
    w.seed_ball(300)
    print(f"   {ind:5s} {NODE_BASIN[ind]:28s} {NODE_DISTANCE_CM[ind]:5.1f} cm -> "
          f"DC transit {w.organs.dc_transit_time():.2f} h")
span = max(NODE_DISTANCE_CM.values()) - min(NODE_DISTANCE_CM.values())
w2 = World(seed=1, indication="BRCA")
w2.seed_ball(300)
w3 = World(seed=1, indication="PDAC")
w3.seed_ball(300)
print(f"   -> a {span:.0f} cm difference changes transit by "
      f"{w2.organs.dc_transit_time() - w3.organs.dc_transit_time():.3f} h, against an "
      f"{PRIMING['dc_transit_base'].value:.0f} h baseline. Distance is NOT rate-limiting.")

print("\n3. Cognate precursor frequency by indication")
for ind in ("CESC", "LUAD", "PDAC"):
    w = World(seed=1, indication=ind)
    w.seed_ball(300)
    kind = "allogeneic cell line" if ind == "CESC" else "autologous neoantigen"
    print(f"   {ind:5s} {w.organs.precursor_frequency:8.1e} ({kind:22s}) -> "
          f"{w.organs.precursors:>12,.0f} precursors in the draining node")

print("\n4. Priming loop timing (LUAD)")
expected = {"dc_arrived_node": None, "effectors_exit_node": None, "effectors_arrive_tumor": None}
w = World(seed=5, indication="LUAD")
w.seed_ball(2000)
w.solve_fields()
w.seed_fibroblasts()
w.seed_immune()
first = {}
for _ in range(5 * 24 * 4):
    w.step()
    ev = [e[1] for e in w.events]
    for k in expected:
        if ev.count(k) and k not in first:
            first[k] = w.t
t_dc = w.organs.dc_transit_time()
print(f"   DC reaches node        {first.get('dc_arrived_node', float('nan')):6.1f} h  "
      f"(transit {t_dc:.1f} h after the first DC leaves)")
print(f"   effectors exit node    {first.get('effectors_exit_node', float('nan')):6.1f} h  "
      f"(+{PRIMING['priming_h'].value:.0f} h priming +{PRIMING['expansion_h'].value:.0f} h expansion)")
print(f"   effectors reach tumour {first.get('effectors_arrive_tumor', float('nan')):6.1f} h  "
      f"(blood circuit ~1 min, so same step)")
ev = [e[1] for e in w.events]
print(f"   -> {ev.count('killed')} tumour cells killed; "
      f"{int((w.tcells.state == EFFECTOR).sum())} effectors in the tumour")

print("\n5. Marrow granulopoiesis (measured: 0.85e9 neutrophils/kg/day, Dancey et al. 1976)")
prod_per_day = MARROW["neutrophil_production"].value * BODY_MASS.value
print(f"   {prod_per_day:.2e} neutrophils/day at {BODY_MASS.value:g} kg")
print(f"   post-mitotic transit {MARROW['postmitotic_transit'].value / 24:.2f} days "
      f"(measured 6.60 +/- 0.03 d)")
print(f"   marrow post-mitotic pool now {w.organs.marrow.count['post-mitotic']:.2e} "
      f"(measured {MARROW['postmitotic_pool'].value * BODY_MASS.value:.2e})")
