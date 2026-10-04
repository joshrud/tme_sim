"""Per-indication tumor configuration and the stochastic somatic-mutation model.

See docs/indications.md for the literature behind every number.

Mutation model (per division, Werner et al. 2020 Nat Commun doi:10.1038/s41467-020-14844-6):
    n_new ~ Poisson(1.14 * tumor_factor)      healthy = 1.14 mut/division; tumors 4-100x
each new mutation independently hits a given driver gene with probability
    driver_target_bp / genome_bp              ~1.5 kb coding / 3.2 Gb = 4.7e-7
and which driver is hit is drawn from the indication's driver weights.

Division-independent mutation is real (Abascal et al. 2021 Nature doi:10.1038/s41586-021-03477-4)
but its magnitude in tumors is unestablished, so `per_hour` defaults to 0.
"""
from dataclasses import dataclass, field

import numpy as np

from .params import P

# Werner et al. 2020: healthy tissue baseline, and the tumor multiplier range
MUT_PER_DIVISION_HEALTHY = P(1.14, "mutations/division", "measured",
                             "healthy haematopoiesis, Werner et al. 2020 doi:10.1038/s41467-020-14844-6")
TUMOR_FACTOR_RANGE = P((4.0, 100.0), "-", "measured", "4-100x healthy in 16 tumors, Werner et al. 2020")
GENOME_BP = P(3.2e9, "bp", "measured", "haploid human genome")
DRIVER_TARGET_BP = P(1500.0, "bp", "assumption",
                     "coding bases per driver gene where a non-synonymous hit is plausible; "
                     "order-of-magnitude, not measured")


@dataclass(frozen=True)
class Indication:
    key: str
    name: str
    line: str                    # reference cell line
    doubling_h: float            # median of all Cellosaurus-reported values
    doubling_range: tuple        # (min, max) reported, to show the spread
    truncal: tuple               # drivers already clonal at t=0
    drivers: dict                # gene -> relative weight of being the next driver hit
    tumor_factor: float          # multiplier on MUT_PER_DIVISION_HEALTHY
    stroma_note: str
    sources: str

    @property
    def mut_per_division(self):
        return MUT_PER_DIVISION_HEALTHY.value * self.tumor_factor


# Driver weights are proportional to reported recurrence; they set *which* gene is hit,
# not how often mutations occur. Frequencies from the TCGA/ICGC papers in docs/indications.md.
INDICATIONS = {
    "PDAC": Indication(
        key="PDAC", name="Pancreatic ductal adenocarcinoma",
        line="PANC-1", doubling_h=29.0, doubling_range=(15.0, 52.0),
        truncal=("KRAS", "TP53", "CDKN2A"),
        drivers={"SMAD4": 20, "RNF43": 8, "ARID1A": 8, "TGFBR2": 6, "GNAS": 5,
                 "RREB1": 5, "PBRM1": 4, "KDM6A": 4},
        tumor_factor=10.0,
        stroma_note="very high: dense desmoplastic stroma, few neoplastic cells",
        sources="Bailey 2016 PMID 26909576; TCGA PDAC 2017 PMID 28810144; Cellosaurus PANC-1"),
    "LUAD": Indication(
        key="LUAD", name="Lung adenocarcinoma",
        line="A549", doubling_h=27.0, doubling_range=(18.0, 40.0),
        truncal=("KRAS", "STK11"),           # A549 genotype; TP53 wild-type in this line
        drivers={"TP53": 46, "KEAP1": 17, "NF1": 11, "RBM10": 8, "MGA": 8,
                 "RIT1": 2, "ERBB2": 3, "MET": 4},
        tumor_factor=40.0,                    # highest TMB of the four (8.9 mut/Mb mean)
        stroma_note="moderate; 52 stromal subtypes reported",
        sources="TCGA LUAD 2014 PMID 25079552; Cellosaurus A549"),
    "OV": Indication(
        key="OV", name="High-grade serous ovarian carcinoma",
        line="Kuramochi", doubling_h=46.0, doubling_range=(26.0, 82.0),
        truncal=("TP53",),                    # 96% of HGSOC tumors
        drivers={"NF1": 12, "BRCA1": 9, "BRCA2": 8, "RB1": 8, "CDK12": 3},
        tumor_factor=15.0,
        stroma_note="moderate; copy-number driven disease",
        sources="TCGA OV 2011 PMID 21720365; Domcke 2013 PMID 23839242 (Kuramochi top-ranked); "
                "Cellosaurus Kuramochi"),
    "BRCA": Indication(
        key="BRCA", name="Breast invasive carcinoma (luminal A reference)",
        line="MCF-7", doubling_h=35.0, doubling_range=(24.0, 80.0),
        truncal=("PIK3CA",),                  # MCF-7 is TP53 wild-type
        drivers={"TP53": 37, "GATA3": 13, "MAP3K1": 9, "CDH1": 8, "MAP2K4": 4,
                 "PTEN": 4, "RB1": 3},
        tumor_factor=8.0,                     # lowest TMB; only 3 genes >10%
        stroma_note="variable by subtype; CAF-S1 rich in TNBC",
        sources="TCGA BRCA 2012 PMID 23000897; Cellosaurus MCF-7"),
    "CESC": Indication(
        key="CESC", name="Cervical carcinoma (HeLa baseline)",
        line="HeLa", doubling_h=20.1, doubling_range=(20.1, 48.0),
        truncal=("HPV18_E6E7",),              # HPV18 integration inactivates TP53/RB pathways
        drivers={"PIK3CA": 14, "TGFBR2": 5, "HLA-A": 5, "EP300": 4, "FBXW7": 4},
        tumor_factor=10.0,
        stroma_note="sparsely characterized",
        sources="TCGA CESC 2017 PMID 28112728; Schwarz 1985 PMID 2983228; "
                "cycle 20.1 h from Puck & Steffen 1963 (Cellosaurus reports 31-48 h; see docs)"),
}

DEFAULT = "CESC"  # keeps existing HeLa runs unchanged


class MutationModel:
    """Per-division somatic mutation accumulation for one indication.

    Tracks, per tumor cell: total mutation count and a bitmask of mutated driver genes.
    Driver fitness effects are small and are an *assumption* - Williams et al. 2016
    (doi:10.1038/ng.3489) found a third of tumors fit a neutral model.
    """

    # assumption: each additional driver shortens the cycle by this fraction, capped
    DRIVER_CYCLE_ADVANTAGE = P(0.03, "-", "assumption",
                               "per-driver cycle shortening; near-neutral default per Williams 2016")
    MAX_ADVANTAGE = P(0.25, "-", "assumption", "cap on total driver advantage")

    def __init__(self, indication, rng, per_hour=0.0):
        self.ind = INDICATIONS[indication] if isinstance(indication, str) else indication
        self.rng = rng
        self.per_hour = per_hour  # division-independent rate; 0 = off (see docs)
        self.genes = list(self.ind.drivers)
        w = np.array([self.ind.drivers[g] for g in self.genes], float)
        self.gene_p = w / w.sum()
        self.p_driver = DRIVER_TARGET_BP.value / GENOME_BP.value
        self.truncal_mask = 0  # truncal drivers live outside the acquired bitmask

    def founder(self, k):
        """Mutation state for k seeded cells: no acquired mutations yet."""
        return np.zeros(k, np.int64), np.zeros(k, np.int64)  # (n_mut, driver_mask)

    def on_division(self, n_mut, driver_mask):
        """Mutations acquired by k daughters at division. Returns updated (n_mut, driver_mask)."""
        k = len(n_mut)
        if k == 0:
            return n_mut, driver_mask
        new = self.rng.poisson(self.ind.mut_per_division, k)
        n_mut = n_mut + new
        # how many of the new mutations land in any driver gene
        n_gene_hits = self.rng.binomial(new, self.p_driver * len(self.genes))
        hit = np.where(n_gene_hits > 0)[0]
        for i in hit:
            for _ in range(int(n_gene_hits[i])):
                g = self.rng.choice(len(self.genes), p=self.gene_p)
                driver_mask[i] |= (1 << g)
        return n_mut, driver_mask

    def n_drivers(self, driver_mask):
        return np.array([bin(int(m)).count("1") for m in np.atleast_1d(driver_mask)])

    def cycle_scale(self, driver_mask):
        """Multiplier on cycle length: more drivers -> slightly faster cycling."""
        adv = np.minimum(self.n_drivers(driver_mask) * self.DRIVER_CYCLE_ADVANTAGE.value,
                         self.MAX_ADVANTAGE.value)
        return 1.0 - adv

    def gene_names(self, driver_mask):
        m = int(driver_mask)
        return [g for i, g in enumerate(self.genes) if m & (1 << i)]


def table():
    """Rows describing every indication, for the parameter table."""
    rows = []
    for k, ind in INDICATIONS.items():
        rows.append(("indication", f"{k}.doubling_h", ind.doubling_h, "h", "measured",
                     f"median of Cellosaurus-reported values for {ind.line} "
                     f"(range {ind.doubling_range[0]:g}-{ind.doubling_range[1]:g} h)"))
        rows.append(("indication", f"{k}.tumor_factor", ind.tumor_factor, "-", "assumption",
                     "placed in the Werner 2020 4-100x range using relative TMB"))
        rows.append(("indication", f"{k}.mut_per_division", round(ind.mut_per_division, 1),
                     "mutations/division", "derived", "1.14 * tumor_factor, Werner et al. 2020"))
    return rows
