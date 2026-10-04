# Four indications: mutation landscape, growth kinetics, and mutation acquisition

Reference material for the `tme_sim` indication layer, covering **PDAC**, **LUAD**, **OV (HGSOC)** and
**BRCA**, plus the **HeLa/CESC** baseline already in the model. Every number is tagged with its evidence
status, matching `tme/params.py`:

| Tag | Meaning |
|---|---|
| **measured** | taken directly from a measurement in the cited work |
| **derived** | computed from measured values |
| **proxy** | measured, but in a different system (other cell line, tissue, or species) |
| **calibrated** | fit so the model reproduces a cited observation |
| **assumption** | no adequate data found; first candidate for fitting or a learned rule |

All references were checked against Europe PMC; PMIDs and DOIs are in [References](#references).

---

## 1. Somatic mutation landscape

### Driver genes by indication

| Indication | Truncal / near-clonal drivers | Other recurrent | Dominant genomic mechanism |
|---|---|---|---|
| **PDAC** | **KRAS**, **TP53**, **CDKN2A**, **SMAD4** | RNF43, ARID1A, TGFBR2, GNAS, RREB1, PBRM1, KDM6A | Point mutation in a fixed 4-gene core, then copy-number loss [1, 2, 3] |
| **LUAD** | **TP53**, **KRAS** *or* **EGFR** (mutually exclusive) | STK11, KEAP1, NF1, RBM10, MGA, RIT1, MET ex14, ERBB2 | High point-mutation load (smoking signature) [4] |
| **OV (HGSOC)** | **TP53** (96% of tumors) | BRCA1, BRCA2, NF1, RB1, CDK12 (each low prevalence) | **Copy-number**, not point mutation: 113 focal CNAs; HR defective in ~50% [5, 6] |
| **BRCA** | **TP53**, **PIK3CA**, **GATA3** (the only genes >10% overall) | MAP3K1, CDH1, MAP2K4, PTEN, RB1; subtype-specific | Strongly subtype-dependent; basal-like resembles HGSOC [7] |
| **CESC** (HeLa) | HPV integration; *PIK3CA*, *TGFBR2*, *HLA-A* | CD274/PDCD1LG2 amplification | Viral oncogene driven (HPV18 E6/E7 in HeLa) [8, 9] |

Two points that matter for modeling:

- **PDAC is the most stereotyped.** Four genes recur in a near-fixed order, so a PDAC cell's starting genotype is
  well defined [1, 3].
- **HGSOC is the least point-mutation driven.** Beyond near-universal *TP53*, recurrent SNVs are rare; the disease
  is shaped by copy-number signatures [5, 6]. A model that only tracks SNVs will under-represent HGSOC evolution.
  Flagged as a known limitation.

### Mutation burden

| Indication | Mutation burden | Status | Source |
|---|---|---|---|
| LUAD | **8.9 mutations/Mb** (mean, TCGA exome) | measured | [4] |
| LUAD | **6.3 mutations/Mb** (median, n=11,855 clinical) | measured | [10] |
| Pan-cancer | **1.5 mutations/Mb** (median non-silent, n=3,083) | measured | [11] |
| Pan-cancer | **3.6 mutations/Mb** (median, n≈100,000) | measured | [10] |
| PDAC, OV, BRCA | low, roughly 1–2/Mb | **assumption** | ordering is well established, but the exact medians are assay-dependent |

> **Caveat.** TMB depends heavily on the assay (panel size, exome vs genome, filtering). Lawrence et al. found
> >1,000-fold variation in median mutation frequency across 27 cancer types, and Chalmers et al. showed TMB rises
> ~2.4-fold between ages 10 and 90 [10, 11]. The model therefore uses the **per-division mutation rate** (§3) as its
> primary parameter and treats TMB only as a cross-check.

---

## 2. In vitro growth and division times

Doubling times were taken from **Cellosaurus**, which aggregates every reported value with its source. Spread
between labs is large (often 2–3 fold), so the table gives the **median of all reported values** and the full range.

| Indication | Reference line | n reports | Median DT | Range | Genomic fidelity to the tumor type |
|---|---|---|---|---|---|
| *(baseline)* | **HeLa** | 3 | **39 h** | 31–48 | CESC; HPV18+ [9] |
| **PDAC** | **PANC-1** | 8 | **29 h** | 15–52 | KRAS/TP53/CDKN2A mutant [12] |
| | MIA PaCa-2 | 6 | 38 h | 26–40 | KRAS/TP53/CDKN2A mutant |
| | BxPC-3 | 3 | 48 h | 48–54 | **KRAS wild-type** (atypical) |
| **LUAD** | **A549** | 11 | **27 h** | 18–40 | KRAS G12S, STK11; TP53 wild-type |
| | NCI-H1975 | 5 | 39 h | 28–42 | EGFR L858R/T790M, TP53 |
| **OV** | **Kuramochi** | 6 | **46 h** | 26–82 | **Top-ranked HGSOC model** [13] |
| | OVCAR-4 | 3 | 41 h | 34–43 | Good HGSOC model [13] |
| | OVCAR-3 | 6 | 51 h | 35–69 | Moderate HGSOC model [13] |
| | *SK-OV-3* | — | — | — | **Poor HGSOC model — avoid** [13] |
| **BRCA** | **MCF-7** | 10 | **35 h** | 24–80 | Luminal A; PIK3CA mut, TP53 wild-type |
| | **MDA-MB-231** | 9 | **31 h** | 25–42 | Basal/TNBC; TP53, KRAS mut |
| | T-47D | 8 | 44 h | 32–60 | Luminal A |
| | SK-BR-3 | 5 | 42 h | 30–60 | HER2-amplified |

**Cell line choice matters more than it looks.** Domcke et al. compared 47 ovarian lines against TCGA HGSOC and
found the two most-cited lines — SK-OV-3 and A2780, together 60% of the literature — are **poorly suited** as HGSOC
models: flat copy-number profiles, no *TP53* mutation, and mutations typical of other subtypes. Kuramochi ranked
top [13]. The model defaults to Kuramochi for OV and records the choice.

> **Note on the HeLa baseline.** The existing model uses a 20.1 h cycle from Puck & Steffen (1963) [14], but
> Cellosaurus reports a 31–48 h median for modern HeLa stocks. The older figure comes from a different lab and era.
> This is an open discrepancy; the current model keeps 20.1 h for continuity, and §5 flags it.

---

## 3. How fast do tumor cells acquire mutations?

**Yes, mutation acquisition is tied to division, and it can be modeled as a stochastic per-division process** —
with one important caveat (below).

### The key measurement

Werner et al. inferred mutation rate **per cell division** from multi-region whole-genome sequencing [15]:

- **1.14 mutations per division** in healthy haematopoiesis
- **1.37 mutations per division** in brain development
- Across 131 biopsies from 16 tumors: **4- to 100-fold increased** mutation rates versus healthy development

So a tumor cell plausibly acquires **~5 to ~115 new mutations per division** genome-wide, with wide variation
between patients. This is the model's primary parameter.

Supporting estimates:
- Adult stem cells accumulate **~40 mutations/year** across small intestine, colon and liver [16].
- Somatic mutation rates are **~2 orders of magnitude above germline** rates [17].

### The caveat: not all mutations need division

Abascal et al. used duplex sequencing to show that **post-mitotic neurons accumulate mutations at a constant rate
through life without dividing**, at rates similar to mitotically active tissues [18]. Differentiated blood and
colon cells carry mutation loads similar to their stem cells despite many more divisions.

**Implication for the model:** a purely division-coupled model is an approximation. The implementation therefore
supports both a per-division term and an optional per-unit-time term, with the time term **off by default** (an
assumption, since its magnitude in tumors is not established).

### Model used

At each division, a daughter acquires

```
n_new ~ Poisson(mu_division)
mu_division = 1.14 * tumor_factor      # tumor_factor in [4, 100], Werner 2020
```

Each new mutation independently hits a driver gene with probability

```
p_driver = driver_target_bp / genome_bp      # ~1.5 kb coding per gene / 3.2 Gb
```

which gives ~4.7x10^-7 per mutation per gene — rare per division, but accumulating across ~10^5 cells over a run.
Which driver is hit is drawn from the indication's driver weights (§1).

**Fitness effects are kept small and are an assumption.** Williams et al. found that 1/3 of tumors across 14 types
fit a *neutral* evolutionary model, with clonal selection occurring before tumor growth began [19]. The default is
therefore near-neutral: drivers give a small cycle-time advantage, flagged as an assumption rather than a
measurement, since per-driver in vivo fitness effects are poorly quantified.

### Per-indication rates used

`tumor_factor` is placed within the Werner range using relative TMB. These placements are **assumptions informed
by measured TMB**, not direct measurements.

| Indication | tumor_factor | mutations/division | Rationale |
|---|---|---|---|
| LUAD | 40 | ~46 | Highest TMB of the four (8.9/Mb mean) [4] |
| OV | 15 | ~17 | Moderate; much of the evolution is copy-number, not captured here [5, 6] |
| PDAC | 10 | ~11 | Low TMB, stereotyped driver set [1, 3] |
| BRCA | 8 | ~9 | Lowest TMB; only 3 genes >10% [7] |
| CESC (HeLa) | 10 | ~11 | Default; HeLa is aneuploid and HPV-driven |

---

## 4. Summary table used by the model

| | PDAC | LUAD | OV (HGSOC) | BRCA |
|---|---|---|---|---|
| Reference line | PANC-1 | A549 | Kuramochi | MCF-7 |
| Doubling time (median) | 29 h | 27 h | 46 h | 35 h |
| Truncal drivers | KRAS, TP53, CDKN2A, SMAD4 | TP53 + KRAS/EGFR | TP53 | TP53/PIK3CA |
| Mutations per division | ~11 | ~46 | ~17 | ~9 |
| Main genomic mechanism | point mutation + CN loss | point mutation | **copy number** | subtype-dependent |
| Stroma | very high (desmoplastic) | moderate | moderate | variable by subtype |

---

## 5. Known gaps and limitations

1. **Copy-number evolution is not modeled.** This matters most for HGSOC, where CNAs dominate [5, 6]. SNV-only
   evolution will understate OV heterogeneity.
2. **Per-driver fitness effects are assumptions.** Defaults are near-neutral, consistent with [19].
3. **Division-independent mutation is off by default**, though it demonstrably occurs [18].
4. **The HeLa cycle length is unresolved**: 20.1 h [14] vs 31–48 h (Cellosaurus). Worth re-deriving from modern
   FUCCI lineage data.
5. **Doubling times are bulk population measurements**, not single-cell cycle times, and conflate cycle length
   with the fraction of cycling cells and cell death. Single-cell FUCCI data would be better for all four lines.
6. **Cell lines are not tumors.** Even the best-matched lines diverge from patient tumors [13], and in vitro
   doubling times exceed in vivo tumor growth rates.

---

## References

1. Bailey P, et al. Genomic analyses identify molecular subtypes of pancreatic cancer. *Nature* 2016. PMID 26909576. doi:10.1038/nature16965
2. Waddell N, et al. Whole genomes redefine the mutational landscape of pancreatic cancer. *Nature* 2015. PMID 25719666. doi:10.1038/nature14169
3. Cancer Genome Atlas Research Network. Integrated genomic characterization of pancreatic ductal adenocarcinoma. *Cancer Cell* 2017. PMID 28810144. doi:10.1016/j.ccell.2017.07.007
4. Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. *Nature* 2014. PMID 25079552. doi:10.1038/nature13385
5. Cancer Genome Atlas Research Network. Integrated genomic analyses of ovarian carcinoma. *Nature* 2011. PMID 21720365. doi:10.1038/nature10166
6. Macintyre G, et al. Copy number signatures and mutational processes in ovarian carcinoma. *Nat Genet* 2018. PMID 30104763. doi:10.1038/s41588-018-0179-8
7. Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature* 2012. PMID 23000897. doi:10.1038/nature11412
8. Cancer Genome Atlas Research Network. Integrated genomic and molecular characterization of cervical cancer. *Nature* 2017. PMID 28112728. doi:10.1038/nature21386
9. Schwarz E, et al. Structure and transcription of human papillomavirus sequences in cervical carcinoma cells. *Nature* 1985. PMID 2983228. doi:10.1038/314111a0
10. Chalmers ZR, et al. Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden. *Genome Med* 2017. PMID 28420421. doi:10.1186/s13073-017-0424-2
11. Lawrence MS, et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. *Nature* 2013. PMID 23770567. doi:10.1038/nature12213
12. Deer EL, et al. Phenotype and genotype of pancreatic cancer cell lines. *Pancreas* 2010. PMID 20418756. doi:10.1097/MPA.0b013e3181c15963
13. Domcke S, et al. Evaluating cell lines as tumour models by comparison of genomic profiles. *Nat Commun* 2013. PMID 23839242. doi:10.1038/ncomms3126
14. Puck TT, Steffen J. Life cycle analysis of mammalian cells I. *Biophys J* 1963. PMID 14062457
15. Werner B, et al. Measuring single cell divisions in human tissues from multi-region sequencing data. *Nat Commun* 2020. PMID 32098957. doi:10.1038/s41467-020-14844-6
16. Blokzijl F, et al. Tissue-specific mutation accumulation in human adult stem cells during life. *Nature* 2016. PMID 27698416. doi:10.1038/nature19768
17. Milholland B, et al. Differences between germline and somatic mutation rates in humans and mice. *Nat Commun* 2017. PMID 28485371. doi:10.1038/ncomms15183
18. Abascal F, et al. Somatic mutation landscapes at single-molecule resolution. *Nature* 2021. PMID 33911282. doi:10.1038/s41586-021-03477-4
19. Williams MJ, et al. Identification of neutral tumor evolution across cancer types. *Nat Genet* 2016. PMID 26780609. doi:10.1038/ng.3489

*Cell line doubling times: Cellosaurus (https://www.cellosaurus.org), which aggregates reported values with their
original sources; medians computed over all reported values per line.*
