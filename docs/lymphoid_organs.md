# Lymphoid organs: bone marrow, thymus, spleen, lymph node

The four off-screen compartments that supply and educate the immune cells entering the tumor, modeled for a
**60-year-old patient**. Evidence tags as in [`indications.md`](indications.md):
**measured / derived / proxy / calibrated / assumption**.

---

## 0. Two corrections worth stating up front

**1. Bone marrow contains HSCs, not hESCs.** Haematopoietic stem cells are tissue-specific adult stem cells
restricted to blood lineages. Human embryonic stem cells are pluripotent cells from the blastocyst and are not
present in an adult. The model implements the HSC → progenitor → monocyte path.

**2. At 60, the thymus is largely involuted, and this changes where naive T cells come from.** TRECs (the DNA
circles left by TCR rearrangement, a direct marker of thymic output) fall **>95% between ages 25 and 60** [1, 2].
More importantly, deuterium labeling, thymectomy and TREC analysis together show that **the adult human naive
T-cell pool is maintained almost exclusively by peripheral T-cell division, not by thymic output** [3] — unlike
mice, where the thymus does sustain it. So for this patient the thymus is a minor source, and the naive
repertoire is maintained by homeostatic proliferation in the periphery. The model reflects that.

---

## 1. Does anatomical distance actually matter?

The user asked to estimate travel times from how close each organ is. Having looked at the numbers, **anatomical
distance is not the rate-limiting step for blood-borne cells**, and the model says so rather than implying a
precision that does not exist:

- **Blood is fast.** The whole blood volume circulates in roughly a minute, so a cell entering circulation reaches
  any organ in ~1 min regardless of whether the distance is 2 cm or 30 cm. Against step sizes of 0.25 h, that is
  effectively instantaneous.
- **Lymph is slow, and it is where the real delay lives.** A dendritic cell carrying tumor antigen to the draining
  node takes on the order of **12–24 h**, dominated by interstitial crawling and entry into lymphatics rather than
  by path length [4, 5].
- **Residence and development times dominate everything.** Neutrophil marrow transit is **6.6 days** [6]; monocytes
  take **1.6 days** post-mitotically before release [7]. These are 2–3 orders of magnitude larger than transit.

So travel time is modeled as **a fixed lymphatic transit plus a small distance term**, and the distance term is
documented as minor. Blood-borne recruitment is treated as same-step.

### Tumor → draining lymph node basin

Distances are the anatomical separation between the primary site and its first-echelon nodal basin, from standard
regional anatomy and nodal staging. They are **approximate**: real distances vary with patient size and tumor
location within the organ.

| Indication | First-echelon nodal basin | Distance | Status |
|---|---|---|---|
| **PDAC** | peripancreatic, celiac, superior mesenteric | ~2 cm | **assumption** (adjacent) |
| **LUAD** | hilar → mediastinal | ~5 cm | **assumption** |
| **CESC** | parametrial → pelvic | ~5 cm | **assumption** |
| **OV** | pelvic → para-aortic | ~8 cm | **assumption** |
| **BRCA** | axillary | ~12 cm | **assumption** |

At a lymphatic flow on the order of mm/s, 12 cm contributes well under an hour — small next to the 12–24 h
baseline transit. **PDAC's node basin is adjacent, yet PDAC is the least T-cell infiltrated of the four**, which
is itself the point: proximity to a lymph node does not produce infiltration. Stromal exclusion dominates
(see [`fibroblasts.md`](fibroblasts.md)).

---

## 2. Organ contents and developmental stages

Total body lymphocytes are on the order of 10¹², and **spleen and lymph nodes — not gut — are the largest
compartments**; only 5–20% of lymphocytes reside in gut, and 1–9% in the lamina propria [8]. Compartment shares
below follow that and the classic distribution surveys [9, 10], scaled to a 70 kg adult.

### Bone marrow — granulopoiesis and monopoiesis

Dancey et al. measured human marrow neutrophil kinetics directly [6]:

| Quantity | Value | Status |
|---|---|---|
| Mitotic pool (promyelocytes + myelocytes) | 2.11 ± 0.36 ×10⁹ /kg | **measured** [6] |
| Post-mitotic pool (metamyelocytes, bands, segs) | 5.59 ± 0.90 ×10⁹ /kg | **measured** [6] |
| Post-mitotic transit time | **6.60 ± 0.03 days** | **measured** [6] |
| Neutrophil production | **0.85 ×10⁹ /kg/day** (~6×10¹⁰/day at 70 kg) | **measured** [6] |
| Monocyte post-mitotic interval before release | **1.6 days** | **measured** [7] |

Whole-body context: total cellular turnover is ~0.33×10¹² cells/day, ~90% of it blood cells [11].

Stages modeled: **HSC → progenitor → mitotic → post-mitotic → released**.

### Thymus — at 60, nearly silent

| Quantity | Value | Status |
|---|---|---|
| Thymic export, young adult (20–25 y) | ~1.6×10⁷ cells/day | **proxy** (modeling literature) |
| TREC decline, 25 → 60 y | **>95%** | **measured** [1, 2] |
| Involution rate | ~3%/year to middle age, ~1%/year after | **measured** [1] |
| Export at 60 y (model) | ~8×10⁵ cells/day | **derived** (young-adult rate × 5%) |
| Naive T-cell maintenance in adults | **peripheral division, not thymic output** | **measured** [3] |

### Spleen and lymph nodes

| Compartment | Lymphocytes | Status | Role in model |
|---|---|---|---|
| **Lymph nodes** (~500–600 in body) | ~40% of total | **assumption** (from [8, 9, 10]) | Priming site: DC presents antigen, naive T cells activate |
| **Spleen** | ~15% of total | **assumption** (from [8, 9, 10]) | Blood filtering, reservoir; secondary priming |
| **Blood** | ~2% of total | **measured** [9, 10] | Transit compartment |

---

## 3. The priming loop, off-screen

This is the mechanism the organs exist to support, and it replaces the current local-only priming with the
anatomically correct route:

```
tumor  --DC carries antigen, CCR7-dependent, 12-24 h-->  draining lymph node
                                                          |
                             naive T cell meets antigen-bearing DC (3-phase priming)
                                                          |
                 effector T cells  <--blood, ~minutes--  proliferate and exit
                                                          |
                                                        tumor
```

- **Antigen transport is DC-mediated and CCR7-dependent.** CD103+/CD141+ DCs traffic tumor antigen to the lymph
  node, where they both directly stimulate CD8 T cells and hand antigen to resident myeloid cells; losing CCR7 in
  these cells causes defective priming and faster tumor growth [5]. Antigen moves inside discrete vesicles and is
  transferred between DC subsets at synapses [12].
- **Priming takes about a day and has three phases:** ~8 h of brief serial DC contacts, then ~12 h of stable
  conjugates with cytokine production, then proliferation on day 2 [13].
- **Clonal expansion** then produces effectors that re-enter blood and reach the tumor within minutes.

Modeled timings:

| Step | Value | Status |
|---|---|---|
| DC tumor → draining node | 18 h + distance/lymph velocity | **proxy** [4, 5] |
| Priming in node (phases 1+2) | 20 h | **measured** [13] |
| Clonal expansion before exit | 48 h | **assumption** (proliferation begins day 2 [13]) |
| Effector node → tumor (blood) | <1 step | **derived** (circulation ~1 min) |
| Expansion factor per primed clone | 100× | **assumption** |

---

## 4. Gaps

1. **Distances are approximate anatomy**, not per-patient measurements, and contribute little anyway.
2. **Organ contents are counts, not agents.** The compartments are deterministic pools with stochastic release;
   there is no spatial structure inside them.
3. **Thymic output at 60 is derived**, by scaling a young-adult rate by the measured TREC decline, rather than
   measured directly at that age.
4. **Expansion factor is an assumption** and strongly affects how many effectors arrive.
5. **No B-cell germinal centre reaction**, no antibody, and no memory compartment.
6. **The spleen is a reservoir only** — it does not yet contribute meaningfully to priming.
7. Per-organ lymphocyte shares are assumptions derived from survey literature rather than from one measurement.

---

## References

1. Palmer DB. The effect of age on thymic function. *Front Immunol* 2013. PMID 24109481. doi:10.3389/fimmu.2013.00316
2. Douek DC, et al. Changes in thymic function with age and during the treatment of HIV infection. *Nature* 1998. PMID 9872319. doi:10.1038/25374
3. den Braber I, et al. Maintenance of peripheral naive T cells is sustained by thymus output in mice but not humans. *Immunity* 2012. PMID 22365666. doi:10.1016/j.immuni.2012.02.006
4. Martín-Fontecha A, et al. Regulation of dendritic cell migration to the draining lymph node: impact on T lymphocyte traffic and priming. *J Exp Med* 2003. PMID 12925677. doi:10.1084/jem.20030448
5. Roberts EW, et al. Critical role for CD103+/CD141+ dendritic cells bearing CCR7 for tumor antigen trafficking and priming of T cell immunity in melanoma. *Cancer Cell* 2016. PMID 27424807. doi:10.1016/j.ccell.2016.06.003
6. Dancey JT, et al. Neutrophil kinetics in man. *J Clin Invest* 1976. PMID 956397. doi:10.1172/JCI108517
7. Patel AA, et al. The fate and lifespan of human monocyte subsets in steady state and systemic inflammation. *J Exp Med* 2017. PMID 28606987. doi:10.1084/jem.20170355
8. Ganusov VV, De Boer RJ. Do most lymphocytes in humans really reside in the gut? *Trends Immunol* 2007. PMID 17964854. doi:10.1016/j.it.2007.08.009
9. Westermann J, Pabst R. Distribution of lymphocyte subsets and natural killer cells in the human body. *Clin Investig* 1992. PMID 1392422. doi:10.1007/BF00184787
10. Blum KS, Pabst R. Lymphocyte numbers and subsets in the human blood. *Immunol Lett* 2007. PMID 17129612. doi:10.1016/j.imlet.2006.10.009
11. Sender R, Milo R. The distribution of cellular turnover in the human body. *Nat Med* 2021. PMID 33432173. doi:10.1038/s41591-020-01182-9
12. Ruhland MK, et al. Visualizing synaptic transfer of tumor antigens among dendritic cells. *Cancer Cell* 2020. PMID 32516589. doi:10.1016/j.ccell.2020.05.002
13. Mempel TR, et al. T-cell priming by dendritic cells in lymph nodes occurs in three distinct phases. *Nature* 2004. PMID 14712275. doi:10.1038/nature02238
14. Pillay J, et al. In vivo labeling with 2H2O reveals a human neutrophil lifespan of 5.4 days. *Blood* 2010. PMID 20410504. doi:10.1182/blood-2010-01-259028
