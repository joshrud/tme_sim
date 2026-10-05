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
| tcell | effector_influx_cap | 40.0 | cells/h | **assumption** | cap on primed effectors entering the tumour per hour, so a single priming burst does not arrive all at once |
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
| fibroblast | diameter | 15.0 | um | **assumption** | modeled as a sphere of roughly tumor-cell volume |
| fibroblast | speed | 10.0 | um/h | **proxy** | human fibroblasts in 3D collagen I, Hakkinen et al. 2011 doi:10.1089/ten.tea.2010.0273; normal fibroblasts, not CAFs |
| fibroblast | division_h | 120.0 | h | **assumption** | CAF turnover in tumors; no direct measurement found |
| fibroblast | myCAF_radius | 50.0 | um | **assumption** | distance to tumor within which a fibroblast is myCAF (TGF-beta proxy) |
| fibroblast | margin | 250.0 | um | **assumption** | fibroblasts stay within this distance of the tumor edge; CAFs are tissue-resident and do not disperse like circulating cells |
| fibroblast | ecm_voxel | 20.0 | um | **numerical** | ECM grid spacing, ~1.3 cell diameters |
| fibroblast | ecm_radius | 30.0 | um | **assumption** | matrix is deposited locally around the cell |
| fibroblast | ecm_deposit | 0.0012 | 1/h | **calibrated** | density added per voxel per hour by a myCAF; set so the PDAC tumor edge reaches 0.72 density at 7 days (>0.7 target; desmoplasia, Whatcott 2015 / Erkan 2008). The same rate for every indication: ordering comes from fibroblast number alone |
| fibroblast | ecm_icaf_scale | 0.3 | - | **assumption** | iCAFs deposit less matrix than myCAFs |
| fibroblast | ecm_decay | 0.002 | 1/h | **assumption** | slow turnover; MMP activity not modeled |
| fibroblast | ecm_barrier | 0.6 | - | **assumption** | T cells cannot move into matrix denser than this. Salmon et al. 2012 report that aligned dense fibers restrict T cells from entering tumor islets; a speed penalty alone reproduces the opposite (slower cells dwell near tumor and kill more) |
| fibroblast | ecm_block | 0.8 | - | **assumption** | fraction of T-cell speed lost in fully dense matrix; the effect is measured (Salmon et al. 2012 doi:10.1172/JCI45817), this functional form is not |
| fibroblast | fraction.PDAC | 0.35 | - | **assumption** | fraction of all cells that are fibroblasts; ordering supported, value is not |
| fibroblast | fraction.OV | 0.15 | - | **assumption** | fraction of all cells that are fibroblasts; ordering supported, value is not |
| fibroblast | fraction.BRCA | 0.15 | - | **assumption** | fraction of all cells that are fibroblasts; ordering supported, value is not |
| fibroblast | fraction.LUAD | 0.1 | - | **assumption** | fraction of all cells that are fibroblasts; ordering supported, value is not |
| fibroblast | fraction.CESC | 0.1 | - | **assumption** | fraction of all cells that are fibroblasts; ordering supported, value is not |
| immune | macrophage.radius | 8.0 | um | **assumption** | macrophages are large; modelled as a sphere |
| immune | macrophage.speed | 2.0 | um/min | **assumption** | slower than T cells, faster than fibroblasts |
| immune | macrophage.lifespan | 480.0 | h | **assumption** | TAMs are long-lived; no clean human tumour figure |
| immune | macrophage.suppress | 0.5 | - | **assumption** | local reduction of CD8 killing; phenomenological |
| immune | dendritic.radius | 7.0 | um | **assumption** | dendritic cell body |
| immune | dendritic.speed | 3.0 | um/min | **assumption** | placeholder |
| immune | dendritic.lifespan | 168.0 | h | **assumption** | placeholder |
| immune | dendritic.suppress | 0.0 | - | **assumption** | local reduction of CD8 killing; phenomenological |
| immune | NK.radius | 5.0 | um | **assumption** | lymphocyte-sized |
| immune | NK.speed | 6.0 | um/min | **assumption** | lymphocyte-like, slower than T cells |
| immune | NK.lifespan | 168.0 | h | **assumption** | placeholder |
| immune | NK.suppress | 0.0 | - | **assumption** | local reduction of CD8 killing; phenomenological |
| immune | B cell.radius | 5.0 | um | **measured** | resting lymphocyte 8-11 um diameter, BioNumbers 108368 |
| immune | B cell.speed | 6.0 | um/min | **proxy** | B-cell motility coefficient ~1/5 of T cells in lymph node, Miller et al. 2002 doi:10.1126/science.1070051 |
| immune | B cell.lifespan | 168.0 | h | **assumption** | placeholder |
| immune | B cell.suppress | 0.0 | - | **assumption** | local reduction of CD8 killing; phenomenological |
| immune | Treg.radius | 4.0 | um | **measured** | resting lymphocyte, BioNumbers 108368 |
| immune | Treg.speed | 8.0 | um/min | **proxy** | T-cell-like motility |
| immune | Treg.lifespan | 336.0 | h | **assumption** | placeholder |
| immune | Treg.suppress | 0.8 | - | **assumption** | local reduction of CD8 killing; phenomenological |
| immune | neutrophil.radius | 5.0 | um | **assumption** | granulocyte |
| immune | neutrophil.speed | 12.0 | um/min | **assumption** | fast-migrating myeloid cell |
| immune | neutrophil.lifespan | 129.60000000000002 | h | **measured** | human neutrophil circulatory lifespan 5.4 days by in vivo 2H2O labelling, Pillay et al. 2010 doi:10.1182/blood-2010-01-259028 - about 10x longer than the <1 day figure from ex vivo labelling that is often quoted |
| immune | neutrophil.suppress | 0.3 | - | **assumption** | local reduction of CD8 killing; phenomenological |
| immune | monocyte.marrow_delay | 38.400000000000006 | h | **measured** | postmitotic interval before classical monocytes leave marrow, Patel et al. 2017 doi:10.1084/jem.20170355 |
| immune | monocyte.circulating | 24.0 | h | **measured** | classical monocytes circulate ~1 day, Patel et al. 2017 |
| immune | monocyte.to_tam | 24.0 | h | **assumption** | time in tissue before becoming a TAM |
| immune | dc_egress | 0.02 | 1/h | **assumption** | fraction of intratumoural DCs that pick up antigen and leave for the draining node each hour; CCR7-dependent trafficking is established (Roberts et al. 2016) but the rate is not measured |
| immune | suppression_radius | 40.0 | um | **assumption** | distance over which a Treg or TAM suppresses CD8 killing |
| immune | PDAC.leukocyte_fraction | 0.12 | - | **assumption** | immune cells per tumour cell; phenotype 'excluded'; reproduces published ordering, not a measurement |
| immune | LUAD.leukocyte_fraction | 0.3 | - | **assumption** | immune cells per tumour cell; phenotype 'inflamed'; reproduces published ordering, not a measurement |
| immune | OV.leukocyte_fraction | 0.2 | - | **assumption** | immune cells per tumour cell; phenotype 'moderate'; reproduces published ordering, not a measurement |
| immune | BRCA.leukocyte_fraction | 0.22 | - | **assumption** | immune cells per tumour cell; phenotype 'variable'; reproduces published ordering, not a measurement |
| immune | CESC.leukocyte_fraction | 0.2 | - | **assumption** | immune cells per tumour cell; phenotype 'moderate'; reproduces published ordering, not a measurement |
| organs | patient_age | 60.0 | years | **assumption** | patient age; set by the user for this model |
| organs | marrow.mitotic_pool | 2110000000.0 | cells/kg | **measured** | promyelocytes + myelocytes, Dancey et al. 1976 doi:10.1172/JCI108517 |
| organs | marrow.postmitotic_pool | 5590000000.0 | cells/kg | **measured** | metamyelocytes, bands, segs, Dancey et al. 1976 |
| organs | marrow.postmitotic_transit | 158.39999999999998 | h | **measured** | 6.60 +/- 0.03 days, Dancey et al. 1976 |
| organs | marrow.neutrophil_production | 850000000.0 | cells/kg/day | **measured** | Dancey et al. 1976 |
| organs | marrow.monocyte_postmitotic | 38.400000000000006 | h | **measured** | classical monocytes leave marrow after 1.6 d, Patel et al. 2017 |
| organs | marrow.hsc_fraction | 0.0001 | - | **assumption** | HSCs as a fraction of marrow cells |
| organs | thymus.export_young | 16000000.0 | cells/day | **proxy** | thymic export at 20-25 y, from the T-cell ageing modelling literature |
| organs | thymus.trec_decline_25_60 | 0.95 | - | **measured** | TRECs fall >95% between 25 and 60 y, Douek 1998 / Palmer 2013 |
| organs | thymus.involution_per_year | 0.03 | 1/year | **measured** | ~3%/year to middle age, ~1%/year after, Palmer 2013 |
| organs | thymus.peripheral_maintenance | True | bool | **measured** | adult human naive T cells are maintained by peripheral division, not thymic output, den Braber et al. 2012 |
| organs | priming.dc_transit_base | 18.0 | h | **proxy** | DC migration to the draining node, CCR7-dependent; Martin-Fontecha 2003, Roberts et al. 2016. Dominated by interstitial crawling, not path length |
| organs | priming.lymph_velocity | 1.0 | mm/s | **proxy** | lymphatic flow; the distance term is minor |
| organs | priming.priming_h | 20.0 | h | **measured** | ~8 h serial DC contacts then ~12 h stable conjugates, Mempel et al. 2004 doi:10.1038/nature02238 |
| organs | priming.expansion_h | 48.0 | h | **assumption** | clonal expansion before exit; proliferation starts day 2 |
| organs | priming.expansion_factor | 100.0 | - | **assumption** | effectors produced per primed naive T cell over ~48 h of division; roughly 6-7 divisions |
| organs | priming.dc_capacity | 10.0 | cells | **assumption** | cognate naive T cells one arriving DC can prime |
| organs | priming.precursor_autologous | 1e-05 | - | **measured** | human naive CD8 precursor frequency per epitope ranges 0.6e-6 to 1.3e-4 across 6 epitopes incl. MART-1 and NY-ESO-1 and is conserved between people; Alanio et al. 2010 Blood doi:10.1182/blood-2009-10-251124. Midpoint used for an autologous tumour neoantigen |
| organs | priming.precursor_allogeneic | 0.07 | - | **proxy** | HeLa is allogeneic to any patient, so its precursor frequency is the alloreactive one (~7%, Suchin et al. 2001), ~4 orders of magnitude above a neoantigen. CESC/HeLa only |
| organs | priming.blood_circuit | 0.016666666666666666 | h | **measured** | whole blood volume circulates in about a minute |
| organs | share.lymph_node | 0.4 | - | **assumption** | from distribution surveys, Westermann & Pabst 1992 |
| organs | share.spleen | 0.15 | - | **assumption** | from distribution surveys, Westermann & Pabst 1992 |
| organs | share.blood | 0.02 | - | **measured** | blood holds ~2% of body lymphocytes, Blum & Pabst 2007 |
| organs | share.marrow | 0.1 | - | **assumption** | remaining compartment estimate |
| organs | thymic_export_at_age | 800000 | cells/day | **derived** | young-adult rate scaled by the measured TREC decline (>95% from 25 to 60) |
| organs | node_distance.PDAC | 2.0 | cm | **assumption** | peripancreatic / celiac basin; distance contributes <1 h vs an 18 h baseline transit |
| organs | node_distance.LUAD | 5.0 | cm | **assumption** | hilar / mediastinal basin; distance contributes <1 h vs an 18 h baseline transit |
| organs | node_distance.CESC | 5.0 | cm | **assumption** | parametrial / pelvic basin; distance contributes <1 h vs an 18 h baseline transit |
| organs | node_distance.OV | 8.0 | cm | **assumption** | pelvic / para-aortic basin; distance contributes <1 h vs an 18 h baseline transit |
| organs | node_distance.BRCA | 12.0 | cm | **assumption** | axillary basin; distance contributes <1 h vs an 18 h baseline transit |
