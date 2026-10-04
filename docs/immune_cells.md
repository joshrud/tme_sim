# Immune cells in the TME

Parameters and evidence for the immune layer beyond the CD8 T cells already in the model
(see [`README`](../README.md#t-cells-v0-naive-and-effector-cd8-t-cells)). Evidence tags are as in
[`indications.md`](indications.md): **measured / derived / proxy / calibrated / assumption**.

---

## 1. The honest starting point

Two things are solidly established and drive the design:

1. **Immune composition differs qualitatively between tumors, in a way that predicts therapy response.**
   Binnewies et al. define tumor immune microenvironment (TIME) subclasses, and parsing them improves prediction of
   immunotherapy responsiveness [1]. The commonly used three-way split is **immune-inflamed**, **immune-excluded**
   (T cells present but held in stroma) and **immune-desert**.
2. **Myeloid cells are usually the dominant immune population**, and macrophage infiltration correlates with poor
   prognosis and chemotherapy resistance across most cancers [2].

What is **not** available at simulation precision: the actual per-indication percentages of each immune cell type.
The pan-cancer atlases ([3] 210 patients/15 cancers for myeloid, [4] 316 donors/21 cancers for T cells, [5] for
lung stroma) report these in figures and supplementary tables rather than in text. **So every composition number
below is an assumption that reproduces the published qualitative ordering, not a measurement.** They are the first
thing to replace.

---

## 2. Cell types modeled

| Type | Color in viewer | Role in model | Key references |
|---|---|---|---|
| **CD8 T cell** | white, red "T" | Kills tumor cells when activated (already implemented) | — |
| **Treg** | white, blue "T" | Suppresses CD8 killing locally | [6] |
| **Macrophage (TAM)** | olive | Long-lived, accumulates in hypoxia, suppresses CD8 | [2, 7, 8] |
| **Dendritic cell** | orange | Rare; licenses naive T cells (enables priming) | [9, 10] |
| **NK cell** | violet | Kills without needing prior antigen exposure | [11] |
| **B cell** | cyan | Clusters into tertiary lymphoid structures at the margin | [12, 13] |
| **Neutrophil** | pale yellow | Short-lived, margin-enriched | [14, 15] |

### Why DCs matter most for what is already built
The existing T-cell layer has naive T cells that **cannot be primed**, because HeLa provides signal 1 without
costimulation, so they go anergic or leave. Dendritic cells are the missing piece: Broz et al. identified a rare
population of activating antigen-presenting cells that is critical for T-cell immunity [9], and expanding
CD103+ DC progenitors at the tumor site improves responses to PD-L1 blockade [10]. Adding DCs closes that loop.

---

## 3. Kinetics

| Quantity | Value | Status | Source |
|---|---|---|---|
| **Classical monocyte: marrow → blood** | 1.6 d postmitotic interval | **measured** | Deuterium labeling in humans [16] |
| **Classical monocyte: time in circulation** | ~1 d | **measured** | [16] |
| **Intermediate / non-classical monocyte lifespan** | ~4 d / ~7 d | **measured** | [16] |
| **Monocyte → TAM differentiation** | 24 h after entering tissue | **assumption** | Monocytes are the main TAM source [7]; the timing is not pinned |
| **TAM lifespan in tissue** | 20 d | **assumption** | TAMs are long-lived; no clean human tumor measurement found |
| **Neutrophil lifespan** | 1 d | **proxy** | Short circulating half-life; tumor neutrophils can live longer [15] |
| **NK / B / DC lifespan** | 7 d | **assumption** | Placeholder |
| **Macrophage speed** | 2 µm/min | **assumption** | Slower than T cells, faster than fibroblasts; not directly sourced |
| **DC speed** | 3 µm/min | **assumption** | — |
| **NK speed** | 6 µm/min | **assumption** | Lymphocyte-like, somewhat slower than T cells |
| **B cell speed** | 6 µm/min | **proxy** | B cells have a motility coefficient ~1/5 that of T cells in lymph node two-photon imaging |
| **Neutrophil speed** | 12 µm/min | **assumption** | Fast migrating myeloid cell |

All immune cells use the same persistent-random-walk motility and the same ECM barrier rule as T cells
(see [`fibroblasts.md`](fibroblasts.md)), so dense stroma excludes them too.

---

## 4. Spatial organization

This is the part the user asked to be driven by empirical data, so each rule cites what it reproduces.

| Cell type | Spatial rule in model | What it reproduces |
|---|---|---|
| **TAM** | Biased toward low oxygen | Macrophages accumulate in hypoxic and necrotic tumor areas, where they adopt distinct pro-angiogenic states [8] |
| **B cell** | Seeded in a few tight clusters at the invasive margin | TLS are ectopic lymphoid aggregates at sites of chronic inflammation, typically peritumoral, and B cells localize within them [12, 13] |
| **Treg** | Co-located with CD8 T cells | CAF-S1 attracts CD4+CD25+ T cells and promotes FOXP3+ Treg differentiation [see `stromal_cells.md`] |
| **Neutrophil** | Margin-enriched | Recruited from circulation; tumor-margin biased [14, 15] |
| **DC** | Sparse, stromal | Rare activating APCs [9] |
| **NK** | Margin-enriched | Circulation-derived |

**Procedural generation.** At setup the model places each population according to these rules and the indication's
composition, giving a different but statistically consistent starting TME on every seed — the "randomized
background individual" idea. The tumor then grows and the immune compartment evolves from there.

---

## 5. Per-indication immune composition

**All assumptions.** They encode the published qualitative picture:

| Indication | Immune phenotype | Leukocyte fraction | CD8 | Treg | TAM | DC | NK | B | Neut |
|---|---|---|---|---|---|---|---|---|---|
| **PDAC** | excluded / desert | 0.12 | 0.08 | 0.07 | **0.50** | 0.02 | 0.03 | 0.10 | 0.20 |
| **LUAD** | inflamed | 0.30 | **0.25** | 0.08 | 0.35 | 0.04 | 0.06 | 0.12 | 0.10 |
| **OV** | moderate | 0.20 | 0.15 | 0.10 | 0.40 | 0.03 | 0.05 | 0.12 | 0.15 |
| **BRCA** | variable by subtype | 0.22 | 0.18 | **0.12** | 0.38 | 0.03 | 0.05 | 0.14 | 0.10 |
| **CESC** | moderate | 0.20 | 0.18 | 0.08 | 0.40 | 0.03 | 0.05 | 0.11 | 0.15 |

Rows after "leukocyte fraction" are fractions *of the immune compartment* and sum to 1.

Rationale for the ordering, which **is** supported:
- **PDAC is myeloid-dominated and T-cell poor**; FAP+ CAFs exclude T cells via CXCL12 [see `stromal_cells.md`].
- **LUAD is the most T-cell infiltrated** of the four, consistent with its high mutational burden
  [see `indications.md`] — neoantigen load associates with CD8 infiltration.
- **Macrophages are the largest immune population in most tumors** [2].
- **Tregs are enriched in TNBC** via CAF-S1 [see `stromal_cells.md`].

---

## 6. Gaps, in priority order

1. **All composition numbers are assumptions.** Replace from the supplementary tables of [3, 4, 5].
2. **Most speeds and lifespans are assumptions**, except monocyte kinetics [16], which are well measured.
3. **No chemokine fields.** Recruitment and positioning use geometric and oxygen-based proxies rather than
   CXCL12/CCL2 gradients. This is the single biggest structural simplification.
4. **TLS are placed, not grown.** Real TLS form through a described neogenesis program [12]; here they are seeded
   as clusters.
5. **Suppression is phenomenological.** Tregs and TAMs reduce CD8 killing through a local multiplier rather than
   through modeled cytokines or checkpoint ligands.
6. **No myeloid-derived suppressor cells**, and no monocyte subsets beyond a single recruited population.

---

## References

1. Binnewies M, et al. Understanding the tumor immune microenvironment (TIME) for effective therapy. *Nat Med* 2018. PMID 29686425. doi:10.1038/s41591-018-0014-x
2. Cassetta L, Pollard JW. Targeting macrophages: therapeutic approaches in cancer. *Nat Rev Drug Discov* 2018. PMID 30361552. doi:10.1038/nrd.2018.169
3. Cheng S, et al. A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells. *Cell* 2021. PMID 33545035. doi:10.1016/j.cell.2021.01.010
4. Zheng L, et al. Pan-cancer single-cell landscape of tumor-infiltrating T cells. *Science* 2021. PMID 34914499. doi:10.1126/science.abe6474
5. Lambrechts D, et al. Phenotype molding of stromal cells in the lung tumor microenvironment. *Nat Med* 2018. PMID 29988129. doi:10.1038/s41591-018-0096-5
6. Togashi Y, et al. Regulatory T cells in cancer immunosuppression — implications for anticancer therapy. *Nat Rev Clin Oncol* 2019. PMID 30705439. doi:10.1038/s41571-019-0175-7
7. Franklin RA, et al. The cellular and molecular origin of tumor-associated macrophages. *Science* 2014. PMID 24812208. doi:10.1126/science.1252510
8. Lewis CE, Pollard JW. Distinct role of macrophages in different tumor microenvironments. *Cancer Res* 2006. PMID 16423985. doi:10.1158/0008-5472.CAN-05-4005
9. Broz ML, et al. Dissecting the tumor myeloid compartment reveals rare activating antigen-presenting cells critical for T cell immunity. *Cancer Cell* 2014. PMID 25446897. doi:10.1016/j.ccell.2014.09.007
10. Salmon H, et al. Expansion and activation of CD103+ dendritic cell progenitors at the tumor site enhances tumor responses to therapeutic PD-L1 and BRAF inhibition. *Immunity* 2016. PMID 27096321. doi:10.1016/j.immuni.2016.03.012
11. Shimasaki N, et al. NK cells for cancer immunotherapy. *Nat Rev Drug Discov* 2020. PMID 31907401. doi:10.1038/s41573-019-0052-1
12. Sautès-Fridman C, et al. Tertiary lymphoid structures in the era of cancer immunotherapy. *Nat Rev Cancer* 2019. PMID 31092904. doi:10.1038/s41568-019-0144-6
13. Helmink BA, et al. B cells and tertiary lymphoid structures promote immunotherapy response. *Nature* 2020. PMID 31942075. doi:10.1038/s41586-019-1922-8
14. Coffelt SB, et al. Neutrophils in cancer: neutral no more. *Nat Rev Cancer* 2016. PMID 27282249. doi:10.1038/nrc.2016.52
15. Jaillon S, et al. Neutrophil diversity and plasticity in tumour progression and therapy. *Nat Rev Cancer* 2020. PMID 32694624. doi:10.1038/s41568-020-0281-y
16. Patel AA, et al. The fate and lifespan of human monocyte subsets in steady state and systemic inflammation. *J Exp Med* 2017. PMID 28606987. doi:10.1084/jem.20170355
17. Miller MJ, et al. Two-photon imaging of lymphocyte motility and antigen response in intact lymph node. *Science* 2002. PMID 12016203. doi:10.1126/science.1070051 *(T-cell motility coefficient 5–6× that of B cells)*
