"""Parameters for the base cell (HeLa), each tagged with its evidence.

status values:
  measured    - taken directly from a HeLa measurement
  derived     - computed from measured values / standard physics
  proxy       - measured, but in a different cell line (needs HeLa data)
  calibrated  - fit so the model reproduces an observed HeLa spheroid feature
  assumption  - no good data found; candidate for fitting / a small learned model
  numerical   - solver setting, not biology
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class P:
    value: float
    unit: str
    status: str
    source: str


# ---------------------------------------------------------------- cell cycle
CYCLE = {
    # Puck & Steffen 1963 Biophys J 3:379, Table II (BioNumbers 106404)
    "G1": P(8.40, "h", "measured", "Puck & Steffen 1963, Biophys J 3:379"),
    "S": P(6.04, "h", "measured", "Puck & Steffen 1963"),
    "G2": P(4.56, "h", "measured", "Puck & Steffen 1963"),
    "M": P(1.10, "h", "measured", "Puck & Steffen 1963"),
    # total = 20.1 h; cf. ~22 h generation time (Posakony et al. 1977 J Cell Biol 74:468),
    # G1 ~6.7 h / S ~8.8 h by live biosensors (Hahn, Jones & Meyer 2009 Cell Cycle 8:1044)
    "cv": P(0.15, "-", "assumption",
            "No verified HeLa-specific intermitotic-time CV found; ~2-3 h SD on ~20 h is "
            "typical of mammalian lines. Fit from FUCCI lineage data when available."),
}
# Mitotic sub-stages (HeLa, H2B-YFP time-lapse): Chakraborty et al. 2008 Dev Cell 15:657,
# Table 1 (BioNumbers 102579-102581). Sum 65 min, consistent with M = 1.1 h above.
# In the model a division event = completion of cytokinesis.
MITOSIS = {
    "nebd_to_metaphase": P(34 / 60, "h", "measured", "34 +/- 6 min, Chakraborty et al. 2008"),
    "metaphase_to_anaphase": P(21 / 60, "h", "measured", "21 +/- 4 min, Chakraborty et al. 2008"),
    "anaphase_to_cytokinesis": P(10 / 60, "h", "measured", "10 +/- 2 min, Chakraborty et al. 2008"),
}
CYCLE_MEAN_H = sum(CYCLE[k].value for k in ("G1", "S", "G2", "M"))  # 20.1 h

# ---------------------------------------------------------------- geometry
VOLUME = {
    # BioNumbers 103719: 1198-4290 um^3, mean 2425; BioNumbers 109386: ~2600 um^3
    "mean": P(2425.0, "um^3", "measured", "BioNumbers 103719 / 109386 (HeLa)"),
    # birth volume = 2/3 of mean so linear growth V_b -> 2 V_b averages ~mean
    "birth": P(1617.0, "um^3", "derived", "2/3 * mean volume (linear growth V_b -> 2V_b)"),
}

# ---------------------------------------------------------------- oxygen
OXYGEN = {
    "D": P(2000.0, "um^2/s", "measured",
           "~water value used for spheroids, Grimes et al. 2014 J R Soc Interface 11:20131124"),
    # (760 - 47 mmHg H2O) * 0.95 (5% CO2) * 0.21 = 142 mmHg; ignores unstirred layer
    "bath": P(142.0, "mmHg", "derived", "21% O2 incubator, humidified, 5% CO2"),
    "solubility": P(1.38e-6, "M/mmHg", "derived", "0.003 mL O2/dL/mmHg (plasma, 37 C)"),
    # 100 pmol/min per 4e4 cells (Seahorse) and 42 pmol/s/1e6 cells (Felser et al. 2014,
    # MiP2014 abstract) both give ~42 amol/cell/s.
    "uptake": P(42e-18, "mol/s/cell", "measured", "Felser et al. 2014 (Oroboros MiP2014), HeLa ROUTINE"),
    "Km": P(1.0, "mmHg", "measured", "'typically below 1 mmHg' - Grimes et al. 2014 ref [40]"),
    "arrest": P(5.0, "mmHg", "assumption", "PhysiCell default proliferation threshold (convention)"),
    "necrosis": P(0.8, "mmHg", "measured", "severe hypoxia <=0.8 mmHg - Grimes et al. 2014 ref [18]"),
    "necrosis_time": P(6.0, "h", "assumption", "mean time to necrosis below threshold; no HeLa data"),
}

# ---------------------------------------------------------------- glucose / lactate
GLUCOSE = {
    "D": P(600.0, "um^2/s", "assumption", "~water value at 37 C; tissue value lower, not HeLa-specific"),
    "bath": P(25.0, "mM", "assumption", "high-glucose DMEM; check the actual culture medium"),
    # 578.8 nmol/1e6 cells/h = 161 amol/cell/s, measured in HEK293 (NOT HeLa)
    "uptake": P(161e-18, "mol/s/cell", "proxy", "Noguchi et al. 2020 Sci Rep doi:10.1038/s41598-020-70000-6, HEK293 Table 2"),
    "Km": P(0.5, "mM", "assumption", "placeholder half-saturation"),
}
LACTATE = {
    "D": P(1000.0, "um^2/s", "assumption", "~small-solute value"),
    "bath": P(0.0, "mM", "assumption", "fresh medium"),
    "per_glucose": P(1.88, "-", "proxy", "lactate/glucose flux ratio, Noguchi et al. 2020 (HEK293)"),
}

# ---------------------------------------------------------------- mitogen (serum/EGF)
# HeLa-Fucci spheroids (~700 um) keep a ~70 um proliferating rim while the interior is
# only mildly hypoxic (HIF-1a+, pimonidazole-) - Onozato et al. 2017 Cancer Sci, doi:10.1111/cas.13178.
# O2 alone cannot explain that rim (core pO2 stays >10 mmHg), so G1 arrest is driven by a
# consumed mitogen; EGF withdrawal makes spheroids quiescent (Laurent et al. 2013 BMC Cancer, doi:10.1186/1471-2407-13-73).
MITOGEN = {
    "D": P(150.0, "um^2/s", "assumption", "EGF-sized protein in water (~1.5e-6 cm^2/s)"),
    "bath": P(1.0, "normalized", "assumption", "medium level = 1"),
    "uptake_rate": P(0.0257, "1/s", "calibrated",
                     "decay length 76 um: spherical solution (R/r)sinh(r/L)/sinh(R/L) = 0.5 "
                     "at 70 um depth in a 700 um spheroid (Onozato 2017 rim)"),
    "arrest": P(0.5, "normalized", "assumption", "G1 arrest below half of medium level"),
}

# ---------------------------------------------------------------- T cells
# Naive T cells need signal 1 (TCR-MHC/peptide) + signal 2 (CD28-CD80/86 costimulation) to
# become effectors (Chen & Flies 2013, Nat Rev Immunol, doi:10.1038/nri3405); naive cells
# have no cytotoxic machinery until they differentiate (Kaech & Cui 2012, doi:10.1038/nri3307).
# Signal 1 without costimulation drives naive T cells into adaptive tolerance / anergy
# (Schwartz 2003, Annu Rev Immunol, doi:10.1146/annurev.immunol.21.120601.141110).
TCELL = {
    "diameter_naive": P(8.0, "um", "measured", "resting lymphocyte 8-11 um, BioNumbers 108368 (low end)"),
    "blast_volume_factor": P(3.0, "-", "measured",
                             "activated T cells enlarge 2-4x in volume (blastogenesis), JoVE 2016 PMC5226628; midpoint"),
    "speed_naive": P(10.3, "um/min", "measured", "median, Mrass et al. 2006 J Exp Med, doi:10.1084/jem.20060710"),
    "turn_angle_median": P(47.5, "deg", "measured",
                           "median turning angle, Mrass 2006 (applied per 1-min substep: assumption)"),
    "speed_effector": P(8.0, "um/min", "measured",
                        "CTL in tumors, 8 +/- 3 and 10 +/- 4 um/min, Boissonnas et al. 2007 J Exp Med, "
                        "doi:10.1084/jem.20061890; 7.9 um/min, Mrass 2006"),
    # Fraction of a polyclonal naive repertoire that recognizes HeLa (allogeneic to any donor).
    "cognate_fraction_naive": P(0.07, "-", "proxy",
                                "~7% alloreactive across an MHC mismatch, Suchin et al. 2001 J Immunol "
                                "(mouse), doi:10.4049/jimmunol.166.2.973"),
    "tumor_costimulation": P(0.0, "bool", "assumption",
                             "HeLa reported B7-1/B7-2 (CD80/86) low - single low-confidence source; tumor "
                             "cells can directly prime naive CD8 T cells in some models (Thompson et al. "
                             "2010 J Exp Med, doi:10.1084/jem.20092454)"),
    "anergy_signal1_h": P(8.0, "h", "assumption",
                          "cumulative cognate contact without costimulation before anergy; analog of the "
                          "~8 h first priming phase (Mempel et al. 2004 Nature, doi:10.1038/nature02238)"),
    # Effector CTL killing (perforin + granzymes)
    "contact_median": P(15.0, "min", "measured",
                        "median CTL-tumor contact, Weigelin et al. 2021 Nat Commun, doi:10.1038/s41467-021-25282-3"),
    "hits_to_kill": P(3.0, "hits", "measured", "~3 serial sublethal hits, Weigelin 2021"),
    "single_hit_lethal": P(0.05, "-", "measured", "~5% of single contacts kill directly, Weigelin 2021"),
    "damage_recovery_median": P(49.0, "min", "measured",
                                "median recovery of sublethal damage, Weigelin 2021 (exponential decay)"),
    "perforin_pore": P(0.5, "min", "measured",
                       "target permeabilized within ~30 s, Lopez et al. 2013 Blood, doi:10.1182/blood-2012-07-446146"),
    "apoptosis_onset": P(2.0, "min", "measured", "caspase-dependent rounding within 2 min, Lopez 2013"),
    "apoptotic_clearance": P(6.0, "h", "assumption",
                             "apoptotic cell persists before removal; no phagocytes in the model yet"),
    "squeeze": P(0.8, "-", "assumption", "T cells may overlap tumor cells by 20% when passing between them"),
    "world_margin": P(200.0, "um", "assumption",
                      "T cells farther than this beyond the tumor edge leave the tissue (exit); naive T "
                      "cells recirculate rather than residing in tissue"),
    "substep": P(1.0, "min", "numerical", "T-cell motility substep"),
    "track_every": P(2.0, "min", "numerical", "T-cell positions recorded for the viewer"),
}
# In vivo context for validation: one CTL kills 2-16 infected cells/day (Halle et al. 2016 Immunity,
# doi:10.1016/j.immuni.2016.01.010); a single tumor-cell kill took ~6 h on average (Breart et al.
# 2008 J Clin Invest, doi:10.1172/jci34388).

# ---------------------------------------------------------------- numerics
NUMERICS = {
    "dt": P(0.25, "h", "numerical", "biology step; fields are quasi-steady (diffusion ~s-min)"),
    # spheres at random close packing (phi=0.64) fill space at 1 cell per V when their
    # spacing is 0.64^(1/3) = 0.86 of the volume-equivalent diameter (cells deform in tissue)
    "packing": P(0.86, "-", "derived", "rest distance = 0.86*(r_i+r_j) -> space-filling tissue"),
    "contact_tol": P(1.15, "-", "numerical", "edge if distance < tol * rest distance"),
    "mech_iters": P(10, "-", "numerical", "overlap-relaxation iterations per step"),
    "repulsion": P(0.25, "-", "numerical", "fraction of overlap removed per iteration"),
    "adhesion": P(0.10, "-", "numerical", "fraction of gap closed per iteration"),
}


def table():
    """All parameters as rows (group, name, value, unit, status, source)."""
    rows = []
    for gname, group in [("cycle", CYCLE), ("mitosis", MITOSIS), ("volume", VOLUME), ("oxygen", OXYGEN),
                         ("glucose", GLUCOSE), ("lactate", LACTATE), ("mitogen", MITOGEN), ("tcell", TCELL),
                         ("numerics", NUMERICS)]:
        for k, p in group.items():
            rows.append((gname, k, p.value, p.unit, p.status, p.source))
    from .fibroblasts import table as fibroblast_table
    from .indications import table as indication_table
    rows += indication_table()
    rows += fibroblast_table()
    return rows
