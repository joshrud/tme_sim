"""Off-screen lymphoid organs: bone marrow, thymus, spleen, draining lymph node.

These compartments are not rendered. They are tracked in memory as **counts by developmental
stage**, because the cells inside them have no meaningful position in the tumour's coordinate
frame. They supply the cells that enter the tumour, and they host the interactions that happen
away from it - most importantly DC antigen presentation and naive T-cell priming.

Modelled for a 60-year-old patient. See docs/lymphoid_organs.md for the evidence.

Two facts drive the design:
  * The adult human naive T-cell pool is maintained by peripheral division, NOT thymic output
    (den Braber et al. 2012, Immunity, doi:10.1016/j.immuni.2012.02.006). TRECs fall >95%
    between 25 and 60, so at 60 the thymus is a minor source.
  * Anatomical distance is not rate-limiting for blood-borne cells: the blood volume
    circulates in about a minute. The real delays are lymphatic transit (12-24 h for a DC
    reaching the draining node) and developmental residence (days).
"""
from collections import OrderedDict

import numpy as np

from .params import P

PATIENT_AGE = P(60.0, "years", "assumption", "patient age; set by the user for this model")
BODY_MASS = P(70.0, "kg", "assumption", "reference adult, for per-kg marrow figures")

# ---------------------------------------------------------------- bone marrow
MARROW = {
    "mitotic_pool": P(2.11e9, "cells/kg", "measured",
                      "promyelocytes + myelocytes, Dancey et al. 1976 doi:10.1172/JCI108517"),
    "postmitotic_pool": P(5.59e9, "cells/kg", "measured",
                          "metamyelocytes, bands, segs, Dancey et al. 1976"),
    "postmitotic_transit": P(6.60 * 24, "h", "measured",
                             "6.60 +/- 0.03 days, Dancey et al. 1976"),
    "neutrophil_production": P(0.85e9, "cells/kg/day", "measured", "Dancey et al. 1976"),
    "monocyte_postmitotic": P(1.6 * 24, "h", "measured",
                              "classical monocytes leave marrow after 1.6 d, Patel et al. 2017"),
    "hsc_fraction": P(1e-4, "-", "assumption", "HSCs as a fraction of marrow cells"),
}

# ---------------------------------------------------------------- thymus
THYMUS = {
    "export_young": P(1.6e7, "cells/day", "proxy",
                      "thymic export at 20-25 y, from the T-cell ageing modelling literature"),
    "trec_decline_25_60": P(0.95, "-", "measured",
                            "TRECs fall >95% between 25 and 60 y, Douek 1998 / Palmer 2013"),
    "involution_per_year": P(0.03, "1/year", "measured",
                             "~3%/year to middle age, ~1%/year after, Palmer 2013"),
    "peripheral_maintenance": P(True, "bool", "measured",
                                "adult human naive T cells are maintained by peripheral division, "
                                "not thymic output, den Braber et al. 2012"),
}

# ---------------------------------------------------------------- secondary lymphoid
LYMPHOCYTES_TOTAL = P(1e12, "cells", "assumption",
                      "order-of-magnitude total body lymphocytes; spleen and lymph nodes are the "
                      "largest compartments, not gut (Ganusov & De Boer 2007)")
COMPARTMENT_SHARE = {           # fraction of all body lymphocytes
    "lymph_node": P(0.40, "-", "assumption", "from distribution surveys, Westermann & Pabst 1992"),
    "spleen": P(0.15, "-", "assumption", "from distribution surveys, Westermann & Pabst 1992"),
    "blood": P(0.02, "-", "measured", "blood holds ~2% of body lymphocytes, Blum & Pabst 2007"),
    "marrow": P(0.10, "-", "assumption", "remaining compartment estimate"),
}

# ---------------------------------------------------------------- priming kinetics
PRIMING = {
    "dc_transit_base": P(18.0, "h", "proxy",
                         "DC migration to the draining node, CCR7-dependent; Martin-Fontecha 2003, "
                         "Roberts et al. 2016. Dominated by interstitial crawling, not path length"),
    "lymph_velocity": P(1.0, "mm/s", "proxy", "lymphatic flow; the distance term is minor"),
    "priming_h": P(20.0, "h", "measured",
                   "~8 h serial DC contacts then ~12 h stable conjugates, Mempel et al. 2004 "
                   "doi:10.1038/nature02238"),
    "expansion_h": P(48.0, "h", "assumption", "clonal expansion before exit; proliferation starts day 2"),
    "expansion_factor": P(100.0, "-", "assumption",
                         "effectors produced per primed naive T cell over ~48 h of division; "
                         "roughly 6-7 divisions"),
    "dc_capacity": P(10.0, "cells", "assumption",
                     "cognate naive T cells one arriving DC can prime"),
    "precursor_autologous": P(1e-5, "-", "measured",
                              "human naive CD8 precursor frequency per epitope ranges 0.6e-6 to "
                              "1.3e-4 across 6 epitopes incl. MART-1 and NY-ESO-1 and is conserved "
                              "between people; Alanio et al. 2010 Blood doi:10.1182/blood-2009-10-251124. "
                              "Midpoint used for an autologous tumour neoantigen"),
    "precursor_allogeneic": P(0.07, "-", "proxy",
                              "HeLa is allogeneic to any patient, so its precursor frequency is the "
                              "alloreactive one (~7%, Suchin et al. 2001), ~4 orders of magnitude "
                              "above a neoantigen. CESC/HeLa only"),
    "blood_circuit": P(1.0 / 60, "h", "measured", "whole blood volume circulates in about a minute"),
}

# Distance from the primary site to its first-echelon nodal basin (cm). Approximate anatomy.
NODE_DISTANCE_CM = {"PDAC": 2.0, "LUAD": 5.0, "CESC": 5.0, "OV": 8.0, "BRCA": 12.0}
NODE_BASIN = {
    "PDAC": "peripancreatic / celiac", "LUAD": "hilar / mediastinal",
    "CESC": "parametrial / pelvic", "OV": "pelvic / para-aortic", "BRCA": "axillary",
}


def thymic_export_per_day(age_years):
    """Thymic output at a given age, scaled from the young-adult rate by the measured TREC decline."""
    young = THYMUS["export_young"].value
    if age_years <= 25:
        return young
    # >95% lost between 25 and 60; interpolate exponentially and continue past 60
    k = -np.log(1 - THYMUS["trec_decline_25_60"].value) / 35.0
    return young * np.exp(-k * (age_years - 25))


class Organ:
    """One off-screen compartment: an ordered set of developmental stages with transit times."""

    def __init__(self, key, name, stages):
        self.key, self.name = key, name
        self.stages = OrderedDict(stages)           # stage -> transit time (h); None = terminal pool
        self.count = OrderedDict((s, 0.0) for s in self.stages)
        self.released = 0.0

    def set(self, stage, n):
        self.count[stage] = float(n)

    def advance(self, hours, inflow=0.0, self_renewing_head=False):
        """Move cells down the stage chain. Returns the number released this step.

        `inflow` enters the first stage that has a transit time. If `self_renewing_head`,
        the first stage is a self-renewing pool (HSCs) held at constant size: it feeds the
        chain at `inflow` without being depleted, rather than accumulating cells forever.
        """
        names = list(self.stages)
        released = 0.0
        for s in reversed(names):
            t = self.stages[s]
            if t is None or t <= 0:
                continue
            moved = self.count[s] * min(hours / t, 1.0)
            self.count[s] -= moved
            i = names.index(s)
            if i + 1 < len(names):
                self.count[names[i + 1]] += moved
            else:
                released += moved
        if inflow:
            if self_renewing_head and len(names) > 1:
                self.count[names[1]] += inflow      # HSC output; HSC pool stays constant
            else:
                first = next((n for n in names if self.stages[n]), names[0])
                self.count[first] += inflow
        self.released += released
        return released

    def total(self):
        return sum(self.count.values())

    def rows(self):
        return [(s, self.count[s]) for s in self.stages]


class LymphoidSystem:
    """The four off-screen organs plus the DC -> lymph node -> effector T cell priming loop."""

    def __init__(self, world, age=None):
        self.w = world
        self.age = PATIENT_AGE.value if age is None else age
        kg = BODY_MASS.value

        self.marrow = Organ("marrow", "Bone marrow", [
            ("HSC", None),                                   # self-renewing pool
            ("progenitor", 48.0),
            ("mitotic", 72.0),                               # promyelocyte / myelocyte
            ("post-mitotic", MARROW["postmitotic_transit"].value),
        ])
        self.marrow.set("HSC", MARROW["hsc_fraction"].value * 1e12)
        self.marrow.set("mitotic", MARROW["mitotic_pool"].value * kg)
        self.marrow.set("post-mitotic", MARROW["postmitotic_pool"].value * kg)
        self.marrow.set("progenitor", MARROW["mitotic_pool"].value * kg * 0.5)

        self.thymus = Organ("thymus", "Thymus", [
            ("DN thymocyte", 168.0), ("DP thymocyte", 168.0), ("SP thymocyte", 120.0),
        ])
        # involuted at 60: pools scaled by the same factor as export
        scale = thymic_export_per_day(self.age) / THYMUS["export_young"].value
        for s, base in (("DN thymocyte", 2e8), ("DP thymocyte", 1e9), ("SP thymocyte", 2e8)):
            self.thymus.set(s, base * scale)

        total_lymph = LYMPHOCYTES_TOTAL.value
        self.spleen = Organ("spleen", "Spleen", [("resident lymphocyte", None)])
        self.spleen.set("resident lymphocyte", total_lymph * COMPARTMENT_SHARE["spleen"].value)

        self.node = Organ("node", "Draining lymph node", [
            ("naive T cell", None), ("antigen-loaded DC", None),
            ("priming", None), ("expanding", None),   # driven by delay queues, not rates
        ])
        # one draining basin's share of nodal lymphocytes (~500 nodes in the body)
        self.node.set("naive T cell", total_lymph * COMPARTMENT_SHARE["lymph_node"].value / 500.0)
        # Only the cognate fraction of the naive repertoire can be primed against this tumour.
        # HeLa is an allogeneic cell line, so it has an artificially large precursor pool; the
        # four real indications are autologous and use the measured per-epitope frequency.
        allo = self.w.indication.key == "CESC"
        self.precursor_frequency = (PRIMING["precursor_allogeneic"].value if allo
                                    else PRIMING["precursor_autologous"].value)
        self.precursors = self.node.count["naive T cell"] * self.precursor_frequency

        self.organs = [self.marrow, self.thymus, self.spleen, self.node]
        self.in_transit = []        # DCs carrying antigen: (arrival_time, n)
        # priming and expansion are fixed delays, not exponential compartments: the 3-phase
        # priming sequence has a characteristic duration, so cells must not exit immediately
        self.priming_q = []         # (ready_time, n)
        self.expanding_q = []       # (ready_time, n)
        self.effectors_ready = 0.0  # primed effectors waiting to enter the tumour
        self.history = []

    # ------------------------------------------------------------------ transit
    def dc_transit_time(self):
        """Tumour -> draining node, for this indication. Distance term is deliberately small."""
        d_cm = NODE_DISTANCE_CM.get(self.w.indication.key, 5.0)
        v_cm_per_h = PRIMING["lymph_velocity"].value * 0.1 * 3600  # mm/s -> cm/h
        return PRIMING["dc_transit_base"].value + d_cm / v_cm_per_h

    def send_dcs(self, n):
        """Antigen-loaded DCs leave the tumour for the draining node."""
        if n > 0:
            self.in_transit.append((self.w.t + self.dc_transit_time(), float(n)))

    # ------------------------------------------------------------------ step
    def step(self, hours):
        w = self.w
        # marrow: steady production, released cells become the pool immune.py recruits from
        prod = MARROW["neutrophil_production"].value * BODY_MASS.value * (hours / 24.0)
        self.marrow.advance(hours, inflow=prod, self_renewing_head=True)
        # thymus: at 60 this is nearly silent
        self.thymus.advance(hours, inflow=thymic_export_per_day(self.age) * (hours / 24.0))
        self.spleen.advance(hours)

        # DCs arriving at the node load antigen onto resident naive T cells
        arrived = sum(n for t_arr, n in self.in_transit if t_arr <= w.t)
        self.in_transit = [(t_arr, n) for t_arr, n in self.in_transit if t_arr > w.t]
        if arrived:
            self.node.count["antigen-loaded DC"] += arrived
            # each DC primes a limited number of *cognate* naive T cells
            recruit = min(arrived * PRIMING["dc_capacity"].value, self.precursors)
            if recruit > 0:
                self.precursors -= recruit
                self.node.count["naive T cell"] -= recruit
                self.priming_q.append((w.t + PRIMING["priming_h"].value, recruit))
            w.events.append((w.t, "dc_arrived_node", int(arrived), -1))

        # fixed-delay queues: priming -> expansion -> exit
        done = sum(n for t_r, n in self.priming_q if t_r <= w.t)
        self.priming_q = [(t_r, n) for t_r, n in self.priming_q if t_r > w.t]
        if done:
            self.expanding_q.append((w.t + PRIMING["expansion_h"].value, done))
        released = sum(n for t_r, n in self.expanding_q if t_r <= w.t)
        self.expanding_q = [(t_r, n) for t_r, n in self.expanding_q if t_r > w.t]
        self.node.count["priming"] = sum(n for _, n in self.priming_q)
        self.node.count["expanding"] = sum(n for _, n in self.expanding_q)
        if released:
            self.effectors_ready += released * PRIMING["expansion_factor"].value
            w.events.append((w.t, "effectors_exit_node", int(released), -1))
        self.history.append((w.t, self.snapshot()))
        return self.effectors_ready

    def take_effectors(self, max_n=None):
        """Effectors leave the node and reach the tumour within one step (blood is fast)."""
        n = self.effectors_ready if max_n is None else min(self.effectors_ready, max_n)
        self.effectors_ready -= n
        return int(n)

    # ------------------------------------------------------------------ reporting
    def snapshot(self):
        return {o.key: {s: v for s, v in o.rows()} for o in self.organs}

    def table(self):
        """Rows for the side table in the viewer: (organ, stage, count)."""
        out = []
        for o in self.organs:
            for s, v in o.rows():
                out.append((o.name, s, v))
        out.append(("Draining lymph node", "cognate precursors", self.precursors))
        out.append(("Draining lymph node", "effectors ready", self.effectors_ready))
        out.append(("In transit", "DCs to node", sum(n for _, n in self.in_transit)))
        return out


def table():
    rows = [("organs", "patient_age", PATIENT_AGE.value, PATIENT_AGE.unit,
             PATIENT_AGE.status, PATIENT_AGE.source)]
    for group, d in (("marrow", MARROW), ("thymus", THYMUS), ("priming", PRIMING)):
        for k, p in d.items():
            rows.append(("organs", f"{group}.{k}", p.value, p.unit, p.status, p.source))
    for k, p in COMPARTMENT_SHARE.items():
        rows.append(("organs", f"share.{k}", p.value, p.unit, p.status, p.source))
    rows.append(("organs", "thymic_export_at_age",
                 round(thymic_export_per_day(PATIENT_AGE.value)), "cells/day", "derived",
                 "young-adult rate scaled by the measured TREC decline (>95% from 25 to 60)"))
    for k, v in NODE_DISTANCE_CM.items():
        rows.append(("organs", f"node_distance.{k}", v, "cm", "assumption",
                     f"{NODE_BASIN[k]} basin; distance contributes <1 h vs an 18 h baseline transit"))
    return rows
