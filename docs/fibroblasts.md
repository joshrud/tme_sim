# Fibroblasts and ECM deposition

Parameters and evidence for the fibroblast layer. The *biology* of CAF states is already covered in
[`stromal_cells.md`](stromal_cells.md); this document covers the **quantities** the simulation needs:
how many fibroblasts, how fast they move and divide, and how much matrix they lay down.

Evidence tags are as in [`indications.md`](indications.md): **measured / derived / proxy / calibrated / assumption**.

---

## 1. How much of a tumor is fibroblast?

This is the weakest-evidenced part of the model, and the table below says so explicitly.

**What is solidly established:**
- Tumors are far from pure. Consensus purity across 21 TCGA types averages **75.3 ± 18.9%** malignant cells, so
  roughly a quarter of a tumor is stroma plus immune on average. LUAD specifically runs **<70% purity** [1].
- **PDAC is the extreme case.** It is defined by desmoplasia: "a paucity of neoplastic cells embedded within a
  dense desmoplastic stroma" [2]. Collagen I and hyaluronan are high in both primary and metastatic lesions, and
  high collagen I predicts shorter survival (median 6.4 vs 14.6 months) [3]. The ratio of αSMA+ area to
  collagen area is an independent prognostic marker across 233 patients [4].
- Fibroblast *states* (myCAF / iCAF / apCAF) and their drivers are well mapped [see `stromal_cells.md`].

**What is not solidly established at the precision a simulation wants:** the actual *fraction* of cells that are
fibroblasts, per indication. The pan-cancer scRNA-seq atlases report this in figures rather than text, and
scRNA-seq fractions are distorted by tissue dissociation, which recovers immune cells far more efficiently than
fibroblasts. So the numbers below are **starting values, not measurements**.

| Indication | Fibroblast fraction used | Status | Basis |
|---|---|---|---|
| **PDAC** | 0.35 | **assumption** | Desmoplasia dominates the tumor mass [2, 3, 4] |
| **OV (HGSOC)** | 0.15 | **assumption** | Moderate stroma; omental metastases are adipocyte-rich [5] |
| **BRCA** | 0.15 | **assumption** | Subtype-dependent; CAF-S1-rich in TNBC [6] |
| **LUAD** | 0.10 | **assumption** | Purity <70%, but much of that is immune, not fibroblast [1, 7] |
| **CESC** (HeLa) | 0.10 | **assumption** | Sparsely characterized |

The *ordering* (PDAC ≫ OV ≈ BRCA > LUAD) is well supported; the absolute values are not. These are the first
parameters to replace with numbers read from the supplementary tables of [7, 8, 9].

---

## 2. Fibroblast kinetics

| Quantity | Value | Status | Source |
|---|---|---|---|
| **Migration speed** | 10 µm/h (0.17 µm/min) | **proxy** | Human fibroblasts migrate in 3D collagen I gels, fastest in cell-derived matrix and collagen, minimally in basement-membrane extract [10]; normal fibroblasts, not CAFs |
| **Speed vs matrix structure** | slower in denser/less aligned matrix | **measured** | Migration speed and invasion distance show a biphasic response to collagen density, best predicted by fiber alignment [11] |
| **Division time** | 120 h (5 days) | **assumption** | Fibroblasts divide in culture every 24–48 h, but CAF turnover in tumors is far slower; no direct measurement found |
| **Diameter** | 15 µm | **assumption** | Fibroblasts are elongated; modeled as a sphere of similar volume to a tumor cell |
| **Lifespan** | not modeled | — | CAFs are long-lived; removal is by the tissue-margin rule only |
| **Tissue margin** | 250 µm beyond the tumor edge | **assumption** | CAFs are tissue-resident and do not disperse like circulating cells; without this they random-walk >1 mm away in 7 days |
| **Self-slowing** | fibroblast speed scaled by local ECM | **derived** | Same matrix-density coupling as for T cells [11, 12]; makes the stroma self-stabilizing |

**For comparison, this makes fibroblasts ~60× slower than T cells** (10.3 µm/min), which is the key reason they
form a stable scaffold rather than circulating through the tumor.

---

## 3. ECM deposition

### The data gap

A direct figure for **collagen secreted per fibroblast per unit time** does not exist in usable units. What is
documented is relative: collagen is a major fibroblast product, >20% of newly synthesized protein in 24 h culture
media; procollagen half-life drops from 120 to 20 minutes when secretion is induced; and activated fibroblasts
can devote over half their protein output to procollagen. Synthesis per cell varies ~8-fold with confluence
alone.

So an absolute rate would be invented precision. **The deposition rate is instead a calibrated parameter.**

### How it is calibrated

ECM is stored on a voxel grid (20 µm voxels) as a normalized density in [0, 1], where 1 means "dense
desmoplastic matrix". The deposition rate is set so that **a PDAC simulation reaches high ECM density
(>0.7) in the stromal compartment within ~7 days**, reproducing the defining feature of PDAC [2, 3]. Other
indications use the same per-fibroblast rate and reach lower densities purely because they have fewer
fibroblasts — so the *rate* is one calibrated number, not five. Measured at 7 days:

| | PDAC | OV | BRCA | LUAD |
|---|---|---|---|---|
| Fibroblasts seeded | 2,790 | 1,272 | 1,272 | 754 |
| ECM density at tumor edge | **0.72** | 0.28 | 0.28 | 0.14 |
| Tumor cells in dense (>0.7) matrix | 72% | 0% | 0% | 0% |
| T-cell speed there | 36% of free | 76% | 76% | 88% |

| Quantity | Value | Status | Basis |
|---|---|---|---|
| ECM deposition per myCAF | 0.0012 density-units/h per voxel | **calibrated** | PDAC tumor edge reaches **0.72** at 7 days (target >0.7) |
| ECM degradation | 0.002 /h | **assumption** | Slow turnover; MMP activity is not modeled |
| Voxel size | 20 µm | **numerical** | ~1.3 cell diameters |
| Deposition radius | 30 µm | **assumption** | Matrix is laid down locally around the cell |

### Why ECM matters downstream

ECM density feeds back on T-cell migration, which is a **measured** coupling, not an invented one. Salmon et al.
imaged T cells in live human lung tumor slices: T cells moved actively in loose fibronectin and collagen regions
but **poorly in dense matrix areas**, with aligned fibers steering them around tumor islands and keeping them out.
Degrading the matrix with collagenase **increased** T-cell contact with cancer cells [12]. A preprint (not
peer-reviewed) reports the same for CTL specifically: motility and killing are both impaired in dense collagen
while the killing machinery itself stays intact [13].

The model implements this in **two** parts:

```
speed    = base_speed * (1 - ecm_block * ecm_density)     ecm_block   = 0.8  (assumption)
passable = ecm_density < ecm_barrier                      ecm_barrier = 0.6  (assumption)
```

so a T cell is slowed in moderate matrix and **cannot enter** matrix denser than 0.6.

> **Why the barrier is necessary — a correction found by validation.** The first implementation had the speed
> penalty only. It produced the *opposite* of the published biology: slowing T cells made them linger near tumor
> cells, so contacts rose and killing nearly doubled. The documented mechanism is not drag but **restriction** —
> Salmon et al. report that aligned dense fibers *restrict T cells from entering tumor islets*, so they accumulate
> in stroma instead [12]. The *form* of both rules is an assumption; the *existence and direction* of the effect
> is measured [12].

**What the barrier does and does not reproduce** (PDAC, 7 days of matrix build-up, then 150 CTLs enter;
3 seeds, immune recruitment disabled so the matrix effect is isolated):

| | T-cell contacts | Kills |
|---|---|---|
| No fibroblasts (control) | 1,856 | 198.7 ± 11.9 |
| With fibroblasts | 1,257 | 203.3 ± 18.2 |
| | **−32%** | **no significant change** |

So the model reproduces **spatial exclusion** — a third fewer tumour contacts — but **not** a reduction in killing
at this CTL dose. That is an honest negative result rather than a success: killing needs ~3 hits inside the
damage-recovery window, and with 150 effectors over 24 h there are still enough contacts to reach that threshold
in the accessible rim.

> An earlier version of this document reported a 9% kill reduction. That came from a **single seed** and did not
> survive replication across three seeds; the effect on kills is within noise. Reducing killing would require
> either a lower effector dose, or matrix that also shields the tumour interior rather than only gating entry.

---

## 4. CAF states in the model

Two states, following the PDAC scheme that generalizes across indications [see `stromal_cells.md`]:

| State | Trigger in model | Behavior | Evidence |
|---|---|---|---|
| **myCAF** | within 50 µm of tumor cells (proxy for TGF-β exposure) | High ECM deposition (1.0×) | myCAFs are αSMA-high, juxtatumoral, TGF-β driven |
| **iCAF** | farther than 50 µm from tumor | Low ECM deposition (0.3×), secretes inflammatory signal | iCAFs are distal and IL-1/JAK-STAT driven |

Using **distance** as the proxy for TGF-β exposure is an assumption. It reproduces the observed spatial
arrangement without requiring a TGF-β field yet; a real TGF-β field is the natural next step, and is noted in
`stromal_cells.md` as a planned addition.

---

## 5. Gaps, in priority order

1. **Fibroblast fractions are assumptions.** Replace with values from the scRNA-seq atlas supplementary tables [7, 8, 9].
2. **ECM deposition rate has no absolute measurement** and is calibrated to a qualitative endpoint.
3. **No TGF-β / IL-1 fields**, so CAF state uses a distance proxy.
4. **CAF division time is an assumption** with no direct source.
5. **Migration speed is a proxy** from non-CAF cells in collagen gels.
6. **ECM is scalar.** Real matrix is anisotropic, and fiber *alignment* predicts motility better than density
   does [11] and steers T cells around tumor islands [12]; only scalar density is modeled.

---

## References

1. Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. *Nat Commun* 2015. PMID 26634437. doi:10.1038/ncomms9971
2. Biffi G, et al. IL1-induced JAK/STAT signaling is antagonized by TGFβ to shape CAF heterogeneity in PDAC. *Cancer Discov* 2019. PMID 30366930. doi:10.1158/2159-8290.CD-18-0710
3. Whatcott CJ, et al. Desmoplasia in primary tumors and metastatic lesions of pancreatic cancer. *Clin Cancer Res* 2015. PMID 25695692. doi:10.1158/1078-0432.CCR-14-1051
4. Erkan M, et al. The activated stroma index is a novel and independent prognostic marker in pancreatic ductal adenocarcinoma. *Clin Gastroenterol Hepatol* 2008. PMID 18639493. doi:10.1016/j.cgh.2008.05.006
5. Nieman KM, et al. Adipocytes promote ovarian cancer metastasis and provide energy for rapid tumor growth. *Nat Med* 2011. PMID 22037646. doi:10.1038/nm.2492
6. Costa A, et al. Fibroblast heterogeneity and immunosuppressive environment in human breast cancer. *Cancer Cell* 2018. PMID 29455927. doi:10.1016/j.ccell.2018.01.011
7. Lambrechts D, et al. Phenotype molding of stromal cells in the lung tumor microenvironment. *Nat Med* 2018. PMID 29988129. doi:10.1038/s41591-018-0096-5
8. Qian J, et al. A pan-cancer blueprint of the heterogeneous tumor microenvironment revealed by single-cell profiling. *Cell Res* 2020. PMID 32561858. doi:10.1038/s41422-020-0355-0
9. Wu SZ, et al. A single-cell and spatially resolved atlas of human breast cancers. *Nat Genet* 2021. PMID 34493872. doi:10.1038/s41588-021-00911-1
10. Hakkinen KM, et al. Direct comparisons of the morphology, migration, cell adhesions, and actin cytoskeleton of fibroblasts in four different three-dimensional extracellular matrices. *Tissue Eng Part A* 2011. PMID 20929283. doi:10.1089/ten.tea.2010.0273
11. Fraley SI, et al. Three-dimensional matrix fiber alignment modulates cell migration and MT1-MMP utility by spatially and temporally directing protrusions. *Sci Rep* 2015. PMID 26423227. doi:10.1038/srep14580
12. Salmon H, et al. Matrix architecture defines the preferential localization and migration of T cells into the stroma of human lung tumors. *J Clin Invest* 2012. PMID 22293174. doi:10.1172/JCI45817
13. Zhao R, et al. Collagen density defines 3D migration of CTLs and their consequent cytotoxicity against tumor cells. *bioRxiv* 2021. doi:10.1101/2021.03.16.435689. **Preprint, not peer-reviewed.**
