| group | parameter | value | unit | status | source |
|---|---|---|---|---|---|
| cycle | G1 | 8.4 | h | **measured** | Puck & Steffen 1963, Biophys J 3:379 |
| cycle | S | 6.04 | h | **measured** | Puck & Steffen 1963 |
| cycle | G2 | 4.56 | h | **measured** | Puck & Steffen 1963 |
| cycle | M | 1.1 | h | **measured** | Puck & Steffen 1963 |
| cycle | cv | 0.15 | - | **assumption** | No verified HeLa-specific intermitotic-time CV found; ~2-3 h SD on ~20 h is typical of mammalian lines. Fit from FUCCI lineage data when available. |
| volume | mean | 2425 | um^3 | **measured** | BioNumbers 103719 / 109386 (HeLa) |
| volume | birth | 1617 | um^3 | **derived** | 2/3 * mean volume (linear growth V_b -> 2V_b) |
| oxygen | D | 2000 | um^2/s | **measured** | ~water value used for spheroids, Grimes et al. 2014 J R Soc Interface 11:20131124 |
| oxygen | bath | 142 | mmHg | **derived** | 21% O2 incubator, humidified, 5% CO2 |
| oxygen | solubility | 1.38e-06 | M/mmHg | **derived** | 0.003 mL O2/dL/mmHg (plasma, 37 C) |
| oxygen | uptake | 4.2e-17 | mol/s/cell | **measured** | Felser et al. 2014 (Oroboros MiP2014), HeLa ROUTINE |
| oxygen | Km | 1 | mmHg | **measured** | 'typically below 1 mmHg' - Grimes et al. 2014 ref [40] |
| oxygen | arrest | 5 | mmHg | **assumption** | PhysiCell default proliferation threshold (convention) |
| oxygen | necrosis | 0.8 | mmHg | **measured** | severe hypoxia <=0.8 mmHg - Grimes et al. 2014 ref [18] |
| oxygen | necrosis_time | 6 | h | **assumption** | mean time to necrosis below threshold; no HeLa data |
| glucose | D | 600 | um^2/s | **assumption** | ~water value at 37 C; tissue value lower, not HeLa-specific |
| glucose | bath | 25 | mM | **assumption** | high-glucose DMEM; check the actual culture medium |
| glucose | uptake | 1.61e-16 | mol/s/cell | **proxy** | Noguchi et al. 2020 Sci Rep doi:10.1038/s41598-020-70000-6, HEK293 Table 2 |
| glucose | Km | 0.5 | mM | **assumption** | placeholder half-saturation |
| lactate | D | 1000 | um^2/s | **assumption** | ~small-solute value |
| lactate | bath | 0 | mM | **assumption** | fresh medium |
| lactate | per_glucose | 1.88 | - | **proxy** | lactate/glucose flux ratio, Noguchi et al. 2020 (HEK293) |
| mitogen | D | 150 | um^2/s | **assumption** | EGF-sized protein in water (~1.5e-6 cm^2/s) |
| mitogen | bath | 1 | normalized | **assumption** | medium level = 1 |
| mitogen | uptake_rate | 0.0257 | 1/s | **calibrated** | decay length 76 um: spherical solution (R/r)sinh(r/L)/sinh(R/L) = 0.5 at 70 um depth in a 700 um spheroid (Onozato 2017 rim) |
| mitogen | arrest | 0.5 | normalized | **assumption** | G1 arrest below half of medium level |
| numerics | dt | 0.25 | h | **numerical** | biology step; fields are quasi-steady (diffusion ~s-min) |
| numerics | packing | 0.86 | - | **derived** | rest distance = 0.86*(r_i+r_j) -> space-filling tissue |
| numerics | contact_tol | 1.15 | - | **numerical** | edge if distance < tol * rest distance |
| numerics | mech_iters | 10 | - | **numerical** | overlap-relaxation iterations per step |
| numerics | repulsion | 0.25 | - | **numerical** | fraction of overlap removed per iteration |
| numerics | adhesion | 0.1 | - | **numerical** | fraction of gap closed per iteration |
