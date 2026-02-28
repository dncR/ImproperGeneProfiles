# ImproperGeneProfiles

Reproducible analysis bundle for the manuscript:

**“Complementary ROC-Derived Indices for Screening Improper Expression Profiles in RNA-Seq Differential Expression Analysis”**

This repository accompanies the revised manuscript and provides the code, processed data, and derived outputs used to:

- run the **simulation study** (one-shot benchmark + sensitivity analysis over selected-set sizes),
- reproduce the **cervical cancer miRNA** case study (DESeq2 inference + ROC-derived screening), and
- regenerate the key **figures** and **saved result objects** used in the paper/Supplementary.

The central framing of the manuscript is reflected here:

- **DESeq2** is treated as the primary *inferential* differential expression (DE) method (p-values / BH-adjusted p-values),
- **cAUC / gAUC / LROC** are treated as *ROC-derived screening/ranking indices* (descriptive; no nominal type-I error control).

## Repository layout

- `R/`
  - `helper_functions.R`: shared helpers (simulation generation, preprocessing, DESeq2 wrappers, ROC-index computation, performance metrics).
  - `simulation_codes.R`: simulation pipeline + figure generation.
  - `real_data_codes.R`: cervical miRNA pipeline + figure generation.
- `data/`
  - `Cervical/`: processed count matrix (`count.csv`) and sample metadata (`condition.csv`) used in the real-data section.
  - `Alzheimer/`: an additional processed dataset kept for completeness (not required for reproducing the revised manuscript figures).
- `saved/`
  - `Simulation/simRes-manuscript.Rda`: saved simulation summary object used for plots/tables.
  - `Cervical/Cervical_results.Rda`: saved real-data results object.
- `figure/`
  - `Simulation/`: simulation figures (TPR/PPV/FPR sensitivity curves; visual diagnostics; one-shot summary).
  - `Cervical/`: cervical figures (scatter diagnostics; volcano; dispersion plot; candidate panels).
- `renv.lock`: pinned R package versions for reproducibility.

## Methods implemented (high level)

### Simulation study

The simulation script generates RNA-seq–like counts under multiple scenarios (sample size and overdispersion), with:

- total features `p = 3000`,
- nominal DE prevalence `π_DEGs` (low-DE and/or high-DE settings, configurable in `simulation_codes.R`),
- **improper genes (IGs)** defined as a subset of true DEGs with non-monotonic “high–low” mixture patterns (see manuscript Methods and comments in code).

Each replicate is processed and evaluated as follows:

1. **Prefiltering** removes low-information features (see `filterCounts()` in `R/helper_functions.R`).
2. **DESeq2 inference** is run on counts; for screening comparisons, features are *ranked* by `padj` (BH-adjusted p-values returned by DESeq2).
3. Counts are normalized and **VST-transformed**, then **ROC-derived indices** (cAUC/gAUC/LROC) are computed per feature and used as *screening scores* (ranked by score).
4. **One-shot benchmark**: the selected-set size is fixed to the nominal design value `K0 = p × π_DEGs` (simulation-only descriptive ranking benchmark).
5. **Sensitivity analysis**: the selected-set size varies over a grid of `K` values to summarize screening trade-offs.

Performance metrics follow the manuscript:

- **TPR** and **PPV** are emphasized in the main text for ranking-based screening behavior across `K`,
- **empirical FPR/specificity** is also computed and provided as a supplementary benchmark (prevalence-dependent behavior can make FPR curves visually overlapped in low-DE settings).

### Cervical cancer miRNA analysis

The real-data script (`R/real_data_codes.R`) reproduces the paper’s analysis flow:

1. load processed `count.csv` and `condition.csv`,
2. prefilter low-information miRNAs,
3. run **DESeq2** (inferential DE),
4. apply **ROC-derived indices** as *screening scores* on VST-transformed values,
5. flag candidate “improper miRNA” profiles using heuristic cutoffs consistent with the revised manuscript narrative:
   - `cAUC < 0.60` and (`gAUC > 0.70` **OR** `LROC > 1.48`) (post-hoc/heuristic screening rule; exploratory).

The DESeq2 “significance” rule reported in the revised manuscript is:

- `|log2FoldChange| > 0.6` and `padj < 0.05` (BH-adjusted).

In code, these thresholds are controlled by variables near the top of `R/real_data_codes.R` (e.g., `logFC_threshold`, `alpha_threshold`) and can be set to match the manuscript exactly.

## How to reproduce (recommended)

### 1) Restore the R environment

From the repository root:

```r
install.packages("renv")
renv::restore()
```

Notes:

- `DESeq2` is a Bioconductor package; `renv::restore()` will install it via Bioconductor as needed.
- The scripts use parallel processing (`parallel::makeCluster()`); you may want to adjust the number of workers for your machine.

### 2) Run simulations

```r
source("R/simulation_codes.R")
```

Outputs:

- figures saved under `figure/Simulation/`
- (optionally) some scripts may also write figures into a manuscript folder path (e.g., `00-manuscript/figure/`) when that folder exists as part of the full paper build; if you cloned this repository alone, you can ignore those extra `ggsave()` targets or create the directory.

Important run-time settings (edit in `R/simulation_codes.R` to match the paper):

- `nSim`: number of replicates per scenario (the revised manuscript reports 1000; you can set this value accordingly),
- `propDE`: DE prevalence (low-DE `0.05`; high-DE `0.30` can be enabled by expanding the scenario grid),
- `nFeat_grid`: selected-set size grid for sensitivity analysis.

### 3) Run the cervical miRNA analysis

```r
source("R/real_data_codes.R")
```

Outputs:

- saved objects under `saved/Cervical/`
- figures under `figure/Cervical/`

## Saved outputs included

To make the repository immediately usable without long reruns, we include:

- `saved/Simulation/simRes-manuscript.Rda`
- `saved/Cervical/Cervical_results.Rda`

These are the objects used to generate the manuscript/Supplementary figures currently stored under `figure/`.

## Citation

If you use this repository, please cite the associated manuscript and DESeq2:

- Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology (2014).

## Contact

For questions about reproducing specific figures/tables, please open an issue on the GitHub repository or contact the corresponding author as listed in the manuscript.
