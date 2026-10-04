| group | parameter | value | unit | status | source |
|---|---|---|---|---|---|
| cycle | G1 | 8.4 | h | **measured** | Puck & Steffen 1963, Biophys J 3:379 |
| cycle | S | 6.04 | h | **measured** | Puck & Steffen 1963 |
| cycle | G2 | 4.56 | h | **measured** | Puck & Steffen 1963 |
| cycle | M | 1.1 | h | **measured** | Puck & Steffen 1963 |
| cycle | cv | 0.15 | - | **assumption** | No verified HeLa-specific intermitotic-time CV found; ~2-3 h SD on ~20 h is typical of mammalian lines. Fit from FUCCI lineage data when available. |
| mitosis | nebd_to_metaphase | 0.5666666666666667 | h | **measured** | 34 +/- 6 min, Chakraborty et al. 2008 |
| mitosis | metaphase_to_anaphase | 0.35 | h | **measured** | 21 +/- 4 min, Chakraborty et al. 2008 |
| mitosis | anaphase_to_cytokinesis | 0.16666666666666666 | h | **measured** | 10 +/- 2 min, Chakraborty et al. 2008 |
| volume | mean | 2425.0 | um^3 | **measured** | BioNumbers 103719 / 109386 (HeLa) |
| volume | birth | 1617.0 | um^3 | **derived** | 2/3 * mean volume (linear growth V_b -> 2V_b) |
| oxygen | D | 2000.0 | um^2/s | **measured** | ~water value used for spheroids, Grimes et al. 2014 J R Soc Interface 11:20131124 |
| oxygen | bath | 142.0 | mmHg | **derived** | 21% O2 incubator, humidified, 5% CO2 |
| oxygen | solubility | 1.38e-06 | M/mmHg | **derived** | 0.003 mL O2/dL/mmHg (plasma, 37 C) |
| oxygen | uptake | 4.2e-17 | mol/s/cell | **measured** | Felser et al. 2014 (Oroboros MiP2014), HeLa ROUTINE |
| oxygen | Km | 1.0 | mmHg | **measured** | 'typically below 1 mmHg' - Grimes et al. 2014 ref [40] |
| oxygen | arrest | 5.0 | mmHg | **assumption** | PhysiCell default proliferation threshold (convention) |
| oxygen | necrosis | 0.8 | mmHg | **measured** | severe hypoxia <=0.8 mmHg - Grimes et al. 2014 ref [18] |
| oxygen | necrosis_time | 6.0 | h | **assumption** | mean time to necrosis below threshold; no HeLa data |
| glucose | D | 600.0 | um^2/s | **assumption** | ~water value at 37 C; tissue value lower, not HeLa-specific |
| glucose | bath | 25.0 | mM | **assumption** | high-glucose DMEM; check the actual culture medium |
| glucose | uptake | 1.61e-16 | mol/s/cell | **proxy** | Noguchi et al. 2020 Sci Rep doi:10.1038/s41598-020-70000-6, HEK293 Table 2 |
| glucose | Km | 0.5 | mM | **assumption** | placeholder half-saturation |
| lactate | D | 1000.0 | um^2/s | **assumption** | ~small-solute value |
| lactate | bath | 0.0 | mM | **assumption** | fresh medium |
| lactate | per_glucose | 1.88 | - | **proxy** | lactate/glucose flux ratio, Noguchi et al. 2020 (HEK293) |
| mitogen | D | 150.0 | um^2/s | **assumption** | EGF-sized protein in water (~1.5e-6 cm^2/s) |
| mitogen | bath | 1.0 | normalized | **assumption** | medium level = 1 |
| mitogen | uptake_rate | 0.0257 | 1/s | **calibrated** | decay length 76 um: spherical solution (R/r)sinh(r/L)/sinh(R/L) = 0.5 at 70 um depth in a 700 um spheroid (Onozato 2017 rim) |
| mitogen | arrest | 0.5 | normalized | **assumption** | G1 arrest below half of medium level |
| tcell | diameter_naive | 8.0 | um | **measured** | resting lymphocyte 8-11 um, BioNumbers 108368 (low end) |
| tcell | blast_volume_factor | 3.0 | - | **measured** | activated T cells enlarge 2-4x in volume (blastogenesis), JoVE 2016 PMC5226628; midpoint |
| tcell | speed_naive | 10.3 | um/min | **measured** | median, Mrass et al. 2006 J Exp Med, doi:10.1084/jem.20060710 |
| tcell | turn_angle_median | 47.5 | deg | **measured** | median turning angle, Mrass 2006 (applied per 1-min substep: assumption) |
| tcell | speed_effector | 8.0 | um/min | **measured** | CTL in tumors, 8 +/- 3 and 10 +/- 4 um/min, Boissonnas et al. 2007 J Exp Med, doi:10.1084/jem.20061890; 7.9 um/min, Mrass 2006 |
| tcell | cognate_fraction_naive | 0.07 | - | **proxy** | ~7% alloreactive across an MHC mismatch, Suchin et al. 2001 J Immunol (mouse), doi:10.4049/jimmunol.166.2.973 |
| tcell | tumor_costimulation | 0.0 | bool | **assumption** | HeLa reported B7-1/B7-2 (CD80/86) low - single low-confidence source; tumor cells can directly prime naive CD8 T cells in some models (Thompson et al. 2010 J Exp Med, doi:10.1084/jem.20092454) |
| tcell | anergy_signal1_h | 8.0 | h | **assumption** | cumulative cognate contact without costimulation before anergy; analog of the ~8 h first priming phase (Mempel et al. 2004 Nature, doi:10.1038/nature02238) |
| tcell | contact_median | 15.0 | min | **measured** | median CTL-tumor contact, Weigelin et al. 2021 Nat Commun, doi:10.1038/s41467-021-25282-3 |
| tcell | hits_to_kill | 3.0 | hits | **measured** | ~3 serial sublethal hits, Weigelin 2021 |
| tcell | single_hit_lethal | 0.05 | - | **measured** | ~5% of single contacts kill directly, Weigelin 2021 |
| tcell | damage_recovery_median | 49.0 | min | **measured** | median recovery of sublethal damage, Weigelin 2021 (exponential decay) |
| tcell | perforin_pore | 0.5 | min | **measured** | target permeabilized within ~30 s, Lopez et al. 2013 Blood, doi:10.1182/blood-2012-07-446146 |
| tcell | apoptosis_onset | 2.0 | min | **measured** | caspase-dependent rounding within 2 min, Lopez 2013 |
| tcell | apoptotic_clearance | 6.0 | h | **assumption** | apoptotic cell persists before removal; no phagocytes in the model yet |
| tcell | squeeze | 0.8 | - | **assumption** | T cells may overlap tumor cells by 20% when passing between them |
| tcell | world_margin | 200.0 | um | **assumption** | T cells farther than this beyond the tumor edge leave the tissue (exit); naive T cells recirculate rather than residing in tissue |
| tcell | substep | 1.0 | min | **numerical** | T-cell motility substep |
| tcell | track_every | 2.0 | min | **numerical** | T-cell positions recorded for the viewer |
| numerics | dt | 0.25 | h | **numerical** | biology step; fields are quasi-steady (diffusion ~s-min) |
| numerics | packing | 0.86 | - | **derived** | rest distance = 0.86*(r_i+r_j) -> space-filling tissue |
| numerics | contact_tol | 1.15 | - | **numerical** | edge if distance < tol * rest distance |
| numerics | mech_iters | 10 | - | **numerical** | overlap-relaxation iterations per step |
| numerics | repulsion | 0.25 | - | **numerical** | fraction of overlap removed per iteration |
| numerics | adhesion | 0.1 | - | **numerical** | fraction of gap closed per iteration |
| indication | PDAC.doubling_h | 29.0 | h | **measured** | median of Cellosaurus-reported values for PANC-1 (range 15-52 h) |
| indication | PDAC.tumor_factor | 10.0 | - | **assumption** | placed in the Werner 2020 4-100x range using relative TMB |
| indication | PDAC.mut_per_division | 11.4 | mutations/division | **derived** | 1.14 * tumor_factor, Werner et al. 2020 |
| indication | LUAD.doubling_h | 27.0 | h | **measured** | median of Cellosaurus-reported values for A549 (range 18-40 h) |
| indication | LUAD.tumor_factor | 40.0 | - | **assumption** | placed in the Werner 2020 4-100x range using relative TMB |
| indication | LUAD.mut_per_division | 45.6 | mutations/division | **derived** | 1.14 * tumor_factor, Werner et al. 2020 |
| indication | OV.doubling_h | 46.0 | h | **measured** | median of Cellosaurus-reported values for Kuramochi (range 26-82 h) |
| indication | OV.tumor_factor | 15.0 | - | **assumption** | placed in the Werner 2020 4-100x range using relative TMB |
| indication | OV.mut_per_division | 17.1 | mutations/division | **derived** | 1.14 * tumor_factor, Werner et al. 2020 |
| indication | BRCA.doubling_h | 35.0 | h | **measured** | median of Cellosaurus-reported values for MCF-7 (range 24-80 h) |
| indication | BRCA.tumor_factor | 8.0 | - | **assumption** | placed in the Werner 2020 4-100x range using relative TMB |
| indication | BRCA.mut_per_division | 9.1 | mutations/division | **derived** | 1.14 * tumor_factor, Werner et al. 2020 |
| indication | CESC.doubling_h | 20.1 | h | **measured** | median of Cellosaurus-reported values for HeLa (range 20.1-48 h) |
| indication | CESC.tumor_factor | 10.0 | - | **assumption** | placed in the Werner 2020 4-100x range using relative TMB |
| indication | CESC.mut_per_division | 11.4 | mutations/division | **derived** | 1.14 * tumor_factor, Werner et al. 2020 |
