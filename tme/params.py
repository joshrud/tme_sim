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
    for gname, group in [("cycle", CYCLE), ("volume", VOLUME), ("oxygen", OXYGEN),
                         ("glucose", GLUCOSE), ("lactate", LACTATE), ("mitogen", MITOGEN),
                         ("numerics", NUMERICS)]:
        for k, p in group.items():
            rows.append((gname, k, p.value, p.unit, p.status, p.source))
    return rows
