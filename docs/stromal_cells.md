# Stromal cell dynamics in the tumor microenvironment

A short primer on the non-immune, non-malignant cells of the TME and how they differ by cancer type. It feeds the
cell-identity layer of `tme_sim`. Every claim is tied to a numbered reference; all references were checked
against Europe PMC (PMID/DOI below).

## 1. What the stroma is
| Component | What it does in tumors | Key refs |
|---|---|---|
| **Cancer-associated fibroblasts (CAFs)** | Most abundant stromal cell in most carcinomas. They deposit and remodel ECM, secrete growth factors and cytokines, and shape immunity. They are heterogeneous and plastic, not one cell type. | [1, 2] |
| **Extracellular matrix (ECM)** | Collagen deposition and lysyl-oxidase crosslinking stiffen tissue. That stiffness boosts integrin/PI3K signaling and invasion. | [3] |
| **Endothelial cells (+ pericytes)** | Tumor vessels are abnormal: leaky, tortuous and poorly perfused. This creates hypoxia and acidity and impairs drug and T-cell delivery. "Normalizing" them can restore function. | [4, 5] |
| **Mesenchymal stem/stromal cells (MSCs)** | Recruited to tumors. In breast models, MSC-derived CCL5 increases cancer-cell motility and metastasis. | [6] |
| **Adipocytes** | Where fat dominates (e.g. omentum), they feed tumor cells fatty acids (FABP4-mediated) and secrete IL-8, promoting homing and growth. | [7] |

## 2. Core CAF dynamics
**Origins.** Mostly from resident fibroblasts or stellate cells (e.g. pancreatic stellate cells) [8]. A minority come
from other lineages, including endothelial-to-mesenchymal transition [9]. Across tissues, fibroblasts arise from
"universal" *Pi16*+ / *Col15a1*+ reservoir states that give rise to activated, disease-associated states [10].

**Activation states are signal-driven and reversible.** This is the key dynamic to model:
- **myCAF** (αSMA-high, ECM-producing) sit right next to tumor cells and are driven by **TGFβ** [8, 11].
- **iCAF** (IL-6, LIF and other inflammatory mediators) sit further from tumor cells and are driven by
  **IL-1 → LIF → JAK/STAT**. TGFβ suppresses this state by down-regulating IL1R1, pushing cells back toward myCAF [11].
- **apCAF** express MHC-II and CD74 but lack classic costimulatory molecules. They can present antigen to
  CD4+ T cells [12].
- Breast cancer has a parallel scheme, CAF-S1 to S4 [13]. The FAP+ CAF-S1 subset splits further, into ECM-myCAF and
  TGFβ-myCAF among others, and these clusters track primary resistance to immunotherapy [14].
- A TGFβ-programmed **LRRC15+** myofibroblast lineage is linked to poor anti-PD-L1 response across six cancer types [15].

**How CAFs shape immunity.**
- **Exclusion:** TGFβ-active stroma keeps CD8+ T cells in the peritumoral stroma. Blocking TGFβ lets them in and
  restores checkpoint response [16, 17]. FAP+ CAF-derived **CXCL12** coats cancer cells and excludes T cells;
  CXCR4 blockade reverses this [18].
- **Suppression:** CAF-S1 recruits CD4+CD25+ T cells (CXCL12), retains them (OX40L, PD-L2, JAM2) and promotes their
  conversion to FOXP3+ Tregs [13]. ECM-myCAF raise PD-1/CTLA-4 on Tregs in a feedback loop [14].
- A pan-cancer TGFβ-associated ECM gene program (C-ECM) predicts PD-1 blockade failure [19].

**CAFs can also restrain tumors.** Deleting αSMA+ myofibroblasts in PDAC models *accelerated* disease, with more Tregs
and reduced survival; anti-CTLA-4 reversed this [20]. Reducing stroma by deleting sonic hedgehog also produced more aggressive tumors [21]. So
"less stroma" is not uniformly good, and a model should keep CAF subtypes separate rather than lump them.

**ECM mechanics.** Collagen crosslinking stiffens the matrix, enhances integrin/PI3K signaling and drives invasion
[3]. This is a mechanical input worth adding to the contact graph later.

## 3. Comparison across indications
| Indication | Stromal burden | Dominant CAF states / stromal features | Effect on T cells / therapy | Refs |
|---|---|---|---|---|
| **Pancreatic (PDAC)** | Very high. Dense desmoplastic stroma, with few neoplastic cells. | myCAF (juxtatumoral, TGFβ), iCAF (distal, IL-1/JAK-STAT), apCAF; LRRC15+ myofibroblasts. | FAP+ CAF CXCL12 excludes T cells. CAF depletion can *worsen* outcome (dual role). | [8, 11, 12, 15, 18, 20, 21, 22] |
| **Breast (esp. TNBC)** | Moderate to high; variable by subtype. | CAF-S1 (FAP+, immunosuppressive) and CAF-S4 (myofibroblastic) accumulate in TNBC. Subclasses separate spatially in scRNA-seq. Collagen crosslinking stiffens tissue. | CAF-S1 builds a Treg-rich niche. ECM-myCAF/TGFβ-myCAF clusters mark primary immunotherapy resistance. | [3, 13, 14, 23] |
| **Colorectal (CRC)** | Variable. High in CMS4 ("mesenchymal", 23% of CRC). | TGFβ-activated stroma, angiogenesis, stromal invasion. Poor-prognosis signatures come mainly from stromal genes. | TGFβ drives T-cell exclusion and blocks Th1 differentiation in MSS CRC metastasis. TGFβ blockade sensitizes to anti-PD-1/PD-L1. | [17, 24, 25] |
| **Lung (NSCLC)** | Moderate. | 52 stromal subtypes, including fibroblasts expressing different collagen sets and endothelial cells that down-regulate immune-cell homing genes. | Stromal subtype markers correlate with survival. Endothelial programs co-vary with checkpoint genes and T-cell activity. | [26] |
| **Head & neck (HNSCC)** | Moderate. | Stromal cells (incl. CAFs) share consistent programs across patients. Malignant p-EMT cells localize to the tumor's leading edge. | Stromal composition refines subtypes. p-EMT predicts nodal metastasis. | [27] |
| **Urothelial (bladder)** | Moderate. Commonly immune-excluded. | TGFβ signaling in fibroblasts. CD8+ T cells trapped in fibroblast- and collagen-rich peritumoral stroma. | Fibroblast TGFβ signature marks atezolizumab non-response. TGFβ + PD-L1 blockade restores infiltration in mice. | [16] |
| **Ovarian (omental metastasis)** | Adipocyte-dominated niche. | Omental adipocytes transfer lipids and secrete IL-8. FABP4 up at the tumor–adipocyte interface. | Drives homing and fuels growth. FABP4 loss impairs metastatic growth. | [7] |
| **Cervical (HeLa's origin)** | Reported, but **sparsely characterized**. | One small scRNA-seq study (3 patients) found 3 CAF subtypes. TCGA found *TGFBR2* and *HLA-A* significantly mutated and *CD274*/*PDCD1LG2* amplified. | Direct CAF–T-cell evidence is limited. The mutations above point at TGFβ and antigen-presentation pathways. | [28, 29] |

**Patterns across indications**
- **TGFβ is the common thread.** TGFβ-programmed myofibroblasts and ECM recur as the stromal correlate of T-cell
  exclusion and checkpoint failure in PDAC, breast, CRC, urothelial and pan-cancer data [14–17, 19].
- **Indications differ mostly in stromal burden and niche.** PDAC is desmoplasia-dominated, CRC is subtype-dependent
  (CMS4), and ovarian omental metastases are adipocyte-dominated.
- **Shared CAF states.** Pan-cancer scRNA-seq finds CAF activation trajectories shared across 10 cancer types, plus
  minor CAFs of endothelial origin [9, 10]. So one CAF state space with indication-specific weights is a reasonable
  modeling choice.

## 4. Implications for `tme_sim`
- **Model CAFs as one identity with a continuum of states, not as fixed types.** A position along the iCAF ↔ myCAF
  axis would be set by local **TGFβ** (toward myCAF) vs **IL-1** (toward iCAF) [11]. apCAF would be a separate branch
  that needs MHC-II. This maps directly onto the planned "identity continuum".
- **New fields:** TGFβ, IL-1/LIF and CXCL12 should become graph fields. CXCL12 and TGFβ act on T-cell motility and
  differentiation [16–18], which matters for the T-cell layer now being built.
- **ECM:** a stiffness/collagen field deposited by myCAFs [3]. It would slow T-cell migration and feed tumor-cell invasion.
- **Per-indication setup:** stromal fraction and CAF-state weights (table above) become the "environment" parameters
  when mapping the model to a cancer type.

## References
1. Sahai E, et al. A framework for advancing our understanding of cancer-associated fibroblasts. *Nat Rev Cancer* 2020. PMID 31980749. doi:10.1038/s41568-019-0238-1
2. Kalluri R. The biology and function of fibroblasts in cancer. *Nat Rev Cancer* 2016. PMID 27550820. doi:10.1038/nrc.2016.73
3. Levental KR, et al. Matrix crosslinking forces tumor progression by enhancing integrin signaling. *Cell* 2009. PMID 19931152. doi:10.1016/j.cell.2009.10.027
4. Jain RK. Normalization of tumor vasculature: an emerging concept in antiangiogenic therapy. *Science* 2005. PMID 15637262. doi:10.1126/science.1104819
5. Carmeliet P, Jain RK. Molecular mechanisms and clinical applications of angiogenesis. *Nature* 2011. PMID 21593862. doi:10.1038/nature10144
6. Karnoub AE, et al. Mesenchymal stem cells within tumour stroma promote breast cancer metastasis. *Nature* 2007. PMID 17914389. doi:10.1038/nature06188
7. Nieman KM, et al. Adipocytes promote ovarian cancer metastasis and provide energy for rapid tumor growth. *Nat Med* 2011. PMID 22037646. doi:10.1038/nm.2492
8. Öhlund D, et al. Distinct populations of inflammatory fibroblasts and myofibroblasts in pancreatic cancer. *J Exp Med* 2017. PMID 28232471. doi:10.1084/jem.20162024
9. Luo H, et al. Pan-cancer single-cell analysis reveals the heterogeneity and plasticity of cancer-associated fibroblasts in the tumor microenvironment. *Nat Commun* 2022. PMID 36333338. doi:10.1038/s41467-022-34395-2
10. Buechler MB, et al. Cross-tissue organization of the fibroblast lineage. *Nature* 2021. PMID 33981032. doi:10.1038/s41586-021-03549-5
11. Biffi G, et al. IL1-induced JAK/STAT signaling is antagonized by TGFβ to shape CAF heterogeneity in pancreatic ductal adenocarcinoma. *Cancer Discov* 2019. PMID 30366930. doi:10.1158/2159-8290.CD-18-0710
12. Elyada E, et al. Cross-species single-cell analysis of pancreatic ductal adenocarcinoma reveals antigen-presenting cancer-associated fibroblasts. *Cancer Discov* 2019. PMID 31197017. doi:10.1158/2159-8290.CD-19-0094
13. Costa A, et al. Fibroblast heterogeneity and immunosuppressive environment in human breast cancer. *Cancer Cell* 2018. PMID 29455927. doi:10.1016/j.ccell.2018.01.011
14. Kieffer Y, et al. Single-cell analysis reveals fibroblast clusters linked to immunotherapy resistance in cancer. *Cancer Discov* 2020. PMID 32434947. doi:10.1158/2159-8290.CD-19-1384
15. Dominguez CX, et al. Single-cell RNA sequencing reveals stromal evolution into LRRC15+ myofibroblasts as a determinant of patient response to cancer immunotherapy. *Cancer Discov* 2020. PMID 31699795. doi:10.1158/2159-8290.CD-19-0644
16. Mariathasan S, et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. *Nature* 2018. PMID 29443960. doi:10.1038/nature25501
17. Tauriello DVF, et al. TGFβ drives immune evasion in genetically reconstituted colon cancer metastasis. *Nature* 2018. PMID 29443964. doi:10.1038/nature25492
18. Feig C, et al. Targeting CXCL12 from FAP-expressing carcinoma-associated fibroblasts synergizes with anti-PD-L1 immunotherapy in pancreatic cancer. *PNAS* 2013. PMID 24277834. doi:10.1073/pnas.1320318110
19. Chakravarthy A, et al. TGF-β-associated extracellular matrix genes link cancer-associated fibroblasts to immune evasion and immunotherapy failure. *Nat Commun* 2018. PMID 30410077. doi:10.1038/s41467-018-06654-8
20. Özdemir BC, et al. Depletion of carcinoma-associated fibroblasts and fibrosis induces immunosuppression and accelerates pancreas cancer with reduced survival. *Cancer Cell* 2014. PMID 24856586. doi:10.1016/j.ccr.2014.04.005
21. Rhim AD, et al. Stromal elements act to restrain, rather than support, pancreatic ductal adenocarcinoma. *Cancer Cell* 2014. PMID 24856585. doi:10.1016/j.ccr.2014.04.021
22. Neesse A, et al. Stromal biology and therapy in pancreatic cancer: ready for clinical translation? *Gut* 2019. PMID 30177543. doi:10.1136/gutjnl-2018-316451
23. Bartoschek M, et al. Spatially and functionally distinct subclasses of breast cancer-associated fibroblasts revealed by single cell RNA sequencing. *Nat Commun* 2018. PMID 30514914. doi:10.1038/s41467-018-07582-3
24. Calon A, et al. Stromal gene expression defines poor-prognosis subtypes in colorectal cancer. *Nat Genet* 2015. PMID 25706628. doi:10.1038/ng.3225
25. Guinney J, et al. The consensus molecular subtypes of colorectal cancer. *Nat Med* 2015. PMID 26457759. doi:10.1038/nm.3967
26. Lambrechts D, et al. Phenotype molding of stromal cells in the lung tumor microenvironment. *Nat Med* 2018. PMID 29988129. doi:10.1038/s41591-018-0096-5
27. Puram SV, et al. Single-cell transcriptomic analysis of primary and metastatic tumor ecosystems in head and neck cancer. *Cell* 2017. PMID 29198524. doi:10.1016/j.cell.2017.10.044
28. Wen S, et al. Analysis of cancer-associated fibroblasts in cervical cancer by single-cell RNA sequencing. *Aging (Albany NY)* 2023. PMID 38157264. doi:10.18632/aging.205353
29. Cancer Genome Atlas Research Network. Integrated genomic and molecular characterization of cervical cancer. *Nature* 2017. PMID 28112728. doi:10.1038/nature21386
