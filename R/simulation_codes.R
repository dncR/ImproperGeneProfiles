# Simulation study
library(tidyverse)
library(parallel)
library(DESeq2)
library(nsROC)
library(patchwork)
library(ggrepel)
library(ggtext)
source("R/helper_functions.R")

# 1 - Simulations ----
## 1.1 Scenarios ----
# We will run 27 combinations in total, each repeated 1000 times.
simulation_scenarios <- expand_grid(
  simID = NA_integer_,
  n = c(100, 300, 500),
  p = 3000,
  propDE = c(.05),
  propIR = .25,
  phi = c(.01, .1, 1),
  sdsignal = 1.5,
  propNZ = .10
) |>
  mutate(
    simID = 1:n()
  )

# Number of replicated simulations for each scenario. (e.g., 250, 1000, etc.). Default is 10 for quick testing.
nSim <- 10

# Register parallel backend via "parallel" package.
cl <- makeCluster(4) # Number of CPUs activated.

# Export required elements to parallel nodes.
clusterExport(
  cl,
  c(
    "createDDSobject", "filterCounts", "diffExp", "selectDEfeatures_DESeq",
    "preProcessCounts", "diffExp_ROC", "selectDEfeatures_ROC", "nLhat",
    "generateCountData", "calculatePerformanceMetrics"
  )
)

# simulation_scenarios <- simulation_scenarios[3, ]
# simulation_scenarios <- simulation_scenarios[1:3, ]
simRes <- list()

## 1.2 Computations ----
for (i in 1:nrow(simulation_scenarios)) {
  scen.i <- simulation_scenarios[i, ]
  printStatus(idx = i, scen.i, nrow(simulation_scenarios))

  # We will force each method to select all DE features (one-shot). "nFeat_oneshot" is the number of DEGs which should be observed in the data theoretically. However, the generated dataset may not have this number of DEGs, and it may slightly deviate from the theoretical value.
  nFeat_oneshot <- ceiling(scen.i$propDE * scen.i$p)

  # Define DE feature cluster sizes, which are flagged as DEGs.
  nFeat_grid <- if (scen.i$propDE == .05) {
    feat_max <- ceiling(2.5 * nFeat_oneshot)
    c(5, 10, seq(from = 25, to = feat_max, by = 25))
  }

  # Perform analysis in parallel cores for data list (DF_List)
  clusterExport(cl, c("nFeat_grid", "nFeat_oneshot", "scen.i"))

  # Set seed within clusters
  clusterSetRNGStream(cl, iseed = 3627)

  # Parallel computing over simulated datasets.
  simRes.i <- parLapply(cl, X = 1:nSim, function(idx = X, n_features = nFeat_grid, n_feat_oneshot = nFeat_oneshot, scen = scen.i, ...) {
    library(tidyverse)
    library(DESeq2)

    # Step 1: Generate RNA-Seq Data and store in a list
    DF <- generateCountData(
      n = scen$n, p = scen$p, K = 2,
      param = 1 / scen$phi, sdsignal = scen$sdsignal,
      DE = scen$propDE, IR = scen$propIR,
      nonzero_prop = 0, min_count = 1,
      allZero.rm = TRUE, tag.samples = TRUE
    )

    # Create DESeqDataSet object.
    dds <- createDDSobject(DF)

    # Step 2: Filter data (Remove low quality genes with near zero variances.
    dds_processed <- filterCounts(dds)

    # We will force each method to select all DE features (one-shot). nFeat is the number of DEGs available in the generated RNA-Sequence dataset.
    n_features <- sort(unique(c(n_features, n_feat_oneshot)))

    # DE Analysis via DESeq2
    dds_diffExp <- diffExp(dds_processed, nonzero = TRUE)

    # Number of valid DEGs in the generated dataset.
    nFeat_valid <- length(dds_diffExp$DE_Genes)

    # Step 3: Transform filtered counts using VST method.
    dds_processed <- preProcessCounts(dds_processed,
      normalize = TRUE, transform = TRUE,
      transformationMethod = "vst", nonzero = TRUE
    )

    # DE Analysis via ROC-based methods.
    dds_diffExp_ROC <- diffExp_ROC(.object = dds_processed)
    # dds_diffExp_ROC <- diffExp_ROC(.object = dds_processed, cluster = cl)

    performanceResults <- lapply(n_features, function(n_feat) {
      calculatePerformanceMetrics(objectDESeq = dds_diffExp, objectROC = dds_diffExp_ROC, .n = n_feat)
    }) |>
      bind_rows()

    performanceResults <- bind_cols(
      tibble(dataID = idx, nGenes_oneshot = n_feat_oneshot, nGenes_valid = nFeat_valid),
      performanceResults
    )
    performanceResults <- bind_cols(scen, performanceResults)

    return(performanceResults)
  })

  simRes[[i]] <- bind_rows(simRes.i)

  # Uncomment the line below to save intermediate results.
  save(simRes, file = "saved/Simulation/simRes-manuscript-rev1.Rda")
}

# 2. Plots ----
# 2.1 Performances ----
load("saved/Simulation/simRes-manuscript-rev1.Rda")

simRes <- lapply(simRes, function(x) {
  x |>
    mutate(
      phi = factor(phi, levels = c(.01, .1, 1), labels = c("Very Low", "Moderate", "Very High")),
      n = factor(n, levels = c(100, 300, 500), labels = c("100", "300", "500")),
      propDE = factor(propDE, levels = c(.05, .30), labels = c("Low", "High")),
      group = factor(group),
      Method = factor(Method, levels = c("DESeq", "AUC", "gAUC", "LROC"))
    ) |>
    mutate(
      FPR = 1 - TNR,
      FNR = 1 - TPR
    )
})

simRes_summary <- lapply(simRes, function(x) {
  x_sub <- x |>
    select(-c(simID:propNZ))

  tmp <- x |>
    select(nGenes_oneshot, nGenes_valid) |>
    mutate(
      nGenes_oneshot = unique(nGenes_oneshot),
      nGenes_valid_min = min(nGenes_valid),
      nGenes_valid_max = max(nGenes_valid)
    ) |>
    select(-nGenes_valid) |>
    slice_head(n = 1)

  x_sub <- x_sub |>
    summarise(
      across(c(TPR, TNR, PPV, NPV, FPR, FNR), mean),
      .by = c(group, Method, nGenes_selected)
    )

  sim_scen <- x |>
    select(c(simID:propNZ)) |>
    slice_head(n = 1)

  bind_cols(sim_scen, tmp, x_sub)
}) |>
  bind_rows()

# True number of DEGs observed in the generated datasets. Values are minimum and maximum number of DEGs over all simulations under each scenario.
shadeArea_valid_DEGs <- simRes_summary |>
  distinct(phi, n, nGenes_valid_min, nGenes_valid_max)

# Overall performances for low DE scenarios.
plotData_overall <- simRes |>
  bind_rows() |>
  filter(nGenes_selected == 150)

fig_low <- ggplot(filter(plotData_overall, propDE == "Low"), aes(x = Method, y = TPR, fill = n)) +
  geom_boxplot(width = .6, outlier.colour = rgb(0, 0, 0, .3)) +
  theme_bw() +
  facet_grid(cols = vars(group), rows = vars(phi)) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(margin = margin(t = 5, b = 5)),
    axis.text.y = element_text(margin = margin(r = 5, l = 5))
  ) +
  guides(fill = guide_legend(title = "Sample size (n) ")) +
  labs(x = "Feature selection strategy", y = "True Positive Rate (TPR)")

print(fig_low)

ggsave(
  filename = "figure/Simulation/simulation_overall_low-oneshot.jpeg", plot = fig_low, device = jpeg,
  width = 20, height = 20, units = "cm", dpi = 320, create.dir = TRUE
)

# Overall true detection rate performances of methods for improper gene profiles under low DE scenarios.
# plotData_IGs <- simRes |>
#   bind_rows() |>
#   filter(nGenes_selected %in% c(10, 25, 50, 150), group == "IGs")

# fig_IGs <- ggplot(filter(plotData_IGs, propDE == "Low"), aes(x = Method, y = TPR, fill = n)) +
#   geom_boxplot(width = .6, outlier.colour = rgb(0, 0, 0, .3)) +
#   theme_bw() +
#   facet_grid(
#     cols = vars(nGenes_selected), rows = vars(phi), scales = "free_y",
#     labeller = labeller(nGenes_selected = function(x) paste0("#Genes = ", x))
#   ) +
#   theme(
#     legend.position = "top",
#     axis.text.x = element_text(margin = margin(t = 5, b = 5)),
#     axis.text.y = element_text(margin = margin(r = 5, l = 5))
#   ) +
#   guides(fill = guide_legend(title = "Sample size (n) ")) +
#   labs(x = "Differential Expression Methods", y = "True Detection Rate (IGs)")

# print(fig_IGs)

# 2.2 Sensitivity analysis ----
# Sensitivity analysis at changing number of features selected.
x_breaks <- c(5, 10, 25, 50, 100, 150, 200, 250, 300, 375)
x_breaks_labels <- c("", "10", "", "50", "100", "150", "200", "250", "300", "375")
x_label_sensitivity <- "Number of selected features (K)"

sensitivity_shaded_area <- list(
  geom_rect(
    data = shadeArea_valid_DEGs,
    aes(xmin = nGenes_valid_min, xmax = nGenes_valid_max, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey70", alpha = .10, colour = NA
  ),
  geom_vline(
    data = shadeArea_valid_DEGs,
    aes(xintercept = nGenes_valid_min),
    inherit.aes = FALSE,
    colour = "grey45", linewidth = .35,
    linetype = "dashed"
  ),
  geom_vline(
    data = shadeArea_valid_DEGs,
    aes(xintercept = nGenes_valid_max),
    inherit.aes = FALSE,
    colour = "grey45", linewidth = .35,
    linetype = "dashed"
  )
)

fig_sensitivity_TPR_DEGs <- ggplot(data = filter(simRes_summary, group == "DEGs"), aes(x = nGenes_selected, y = TPR, colour = Method, linetype = Method)) +
  sensitivity_shaded_area +
  geom_point() +
  geom_line(aes(group = Method)) +
  theme_bw(base_size = 12) +
  facet_grid(
    cols = vars(n), rows = vars(phi),
    labeller = labeller(n = function(x) paste0("n = ", x))
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5, b = 5)),
    axis.text.y = element_text(margin = margin(r = 5, l = 5))
  ) +
  guides(fill = guide_legend(title = "Sample size (n) ")) +
  labs(x = x_label_sensitivity, y = "True Positive Rate (TPR) - DEGs") +
  scale_x_continuous(breaks = x_breaks, labels = x_breaks_labels)

ggsave(
  filename = "figure/Simulation/simulation_sensitivity_TPR_DEGs.jpeg", plot = fig_sensitivity_TPR_DEGs, device = jpeg,
  width = 24, height = 20, units = "cm", dpi = 320, create.dir = TRUE
)

fig_sensitivity_TPR_IGs <- ggplot(data = filter(simRes_summary, group == "IGs"), aes(x = nGenes_selected, y = TPR, colour = Method, linetype = Method)) +
  sensitivity_shaded_area +
  geom_point() +
  geom_line(aes(group = Method)) +
  theme_bw(base_size = 12) +
  facet_grid(
    cols = vars(n), rows = vars(phi),
    labeller = labeller(n = function(x) paste0("n = ", x))
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5, b = 5)),
    axis.text.y = element_text(margin = margin(r = 5, l = 5))
  ) +
  guides(fill = guide_legend(title = "Sample size (n) ")) +
  labs(x = x_label_sensitivity, y = "True Positive Rate (TPR) - IGs") +
  scale_x_continuous(breaks = x_breaks, labels = x_breaks_labels) +
  lims(y = c(0, 1))

ggsave(
  filename = "figure/Simulation/simulation_sensitivity_TPR_IGs.jpeg", plot = fig_sensitivity_TPR_IGs, device = jpeg,
  width = 24, height = 20, units = "cm", dpi = 320, create.dir = TRUE
)

fig_sensitivity_PPV_DEGs <- ggplot(data = filter(simRes_summary, group == "DEGs"), aes(x = nGenes_selected, y = PPV, colour = Method, linetype = Method)) +
  sensitivity_shaded_area +
  geom_point() +
  geom_line(aes(group = Method)) +
  theme_bw(base_size = 12) +
  facet_grid(
    cols = vars(n), rows = vars(phi),
    labeller = labeller(n = function(x) paste0("n = ", x))
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5, b = 5)),
    axis.text.y = element_text(margin = margin(r = 5, l = 5))
  ) +
  guides(fill = guide_legend(title = "Sample size (n) ")) +
  labs(x = x_label_sensitivity, y = "Positive Predictive Value (PPV) - DEGs") +
  scale_x_continuous(breaks = x_breaks, labels = x_breaks_labels) +
  lims(y = c(0, 1))

ggsave(
  filename = "figure/Simulation/simulation_sensitivity_PPV_DEGs.jpeg", plot = fig_sensitivity_PPV_DEGs, device = jpeg,
  width = 24, height = 20, units = "cm", dpi = 320, create.dir = TRUE
)

fig_sensitivity_PPV_IGs <- ggplot(data = filter(simRes_summary, group == "IGs"), aes(x = nGenes_selected, y = PPV, colour = Method, linetype = Method)) +
  sensitivity_shaded_area +
  geom_point() +
  geom_line(aes(group = Method)) +
  theme_bw(base_size = 12) +
  facet_grid(
    cols = vars(n), rows = vars(phi),
    labeller = labeller(n = function(x) paste0("n = ", x))
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5, b = 5)),
    axis.text.y = element_text(margin = margin(r = 5, l = 5))
  ) +
  guides(fill = guide_legend(title = "Sample size (n) ")) +
  labs(x = x_label_sensitivity, y = "Positive Predictive Value (PPV) - IGs") +
  scale_x_continuous(breaks = x_breaks, labels = x_breaks_labels) +
  lims(y = c(0, .3))

ggsave(
  filename = "figure/Simulation/simulation_sensitivity_PPV_IGs.jpeg", plot = fig_sensitivity_PPV_IGs, device = jpeg,
  width = 24, height = 20, units = "cm", dpi = 320, create.dir = TRUE
)

fig_sensitivity_FPR_DEGs <- ggplot(data = filter(simRes_summary, group == "DEGs"), aes(x = nGenes_selected, y = FPR, colour = Method, linetype = Method)) +
  sensitivity_shaded_area +
  geom_point() +
  geom_line(aes(group = Method)) +
  theme_bw(base_size = 12) +
  facet_grid(
    cols = vars(n), rows = vars(phi),
    labeller = labeller(n = function(x) paste0("n = ", x))
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5, b = 5)),
    axis.text.y = element_text(margin = margin(r = 5, l = 5))
  ) +
  guides(fill = guide_legend(title = "Sample size (n) ")) +
  labs(x = x_label_sensitivity, y = "False Positive Rate (FPR) - DEGs") +
  scale_x_continuous(breaks = x_breaks, labels = x_breaks_labels) +
  lims(y = c(0, .15))

ggsave(
  filename = "figure/Simulation/simulation_sensitivity_FPR_DEGs.jpeg", plot = fig_sensitivity_FPR_DEGs, device = jpeg,
  width = 24, height = 20, units = "cm", dpi = 320, create.dir = TRUE
)

fig_sensitivity_FPR_IGs <- ggplot(data = filter(simRes_summary, group == "IGs"), aes(x = nGenes_selected, y = FPR, colour = Method, linetype = Method)) +
  sensitivity_shaded_area +
  geom_point() +
  geom_line(aes(group = Method)) +
  theme_bw(base_size = 12) +
  facet_grid(
    cols = vars(n), rows = vars(phi),
    labeller = labeller(n = function(x) paste0("n = ", x))
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5, b = 5)),
    axis.text.y = element_text(margin = margin(r = 5, l = 5))
  ) +
  guides(fill = guide_legend(title = "Sample size (n) ")) +
  labs(x = x_label_sensitivity, y = "False Positive Rate (FPR) - IGs") +
  scale_x_continuous(breaks = x_breaks, labels = x_breaks_labels) +
  lims(y = c(0, .2))

ggsave(
  filename = "figure/Simulation/simulation_sensitivity_FPR_IGs.jpeg", plot = fig_sensitivity_FPR_IGs, device = jpeg,
  width = 24, height = 20, units = "cm", dpi = 320, create.dir = TRUE
)

# Overall performances for high DE scenarios.
# fig_high <- ggplot(filter(simRes_long, propDE != "Low"), aes(x = Method, y = Value, fill = n)) +
#   geom_boxplot(width = .6, outlier.colour = rgb(0, 0, 0, .3)) +
#   theme_bw() +
#   facet_grid(cols = vars(group), rows = vars(phi), scales = "free_y") +
#   theme(
#     legend.position = "top",
#     axis.text.x = element_text(margin = margin(t = 5, b = 5)),
#     axis.text.y = element_text(margin = margin(r = 5, l = 5))
#   ) +
#   guides(fill = guide_legend(title = "Sample size (n) ")) +
#   labs(x = "Differential Expression Methods", y = "True Detection Rate")

# Save figure
# ggsave(filename = "figure/Simulation/simRes_overall_high.png", plot = fig_high, device = png,
#        width = 18, height = 22, units = "cm", dpi = 320, create.dir = TRUE)

# ggsave(filename = "document/manuscript/figure/simRes_overall_high.png", plot = fig_high, device = png,
#        width = 18, height = 22, units = "cm", dpi = 320, create.dir = TRUE)

# ggsave(filename = "figure/Simulation/simRes_overall.svg", plot = fig, device = svg,
#        width = 18, height = 26, units = "cm", dpi = 320, create.dir = TRUE)

# 2.3 IG Inspection ----
# Step 1: Generate RNA-Seq Data and store in a list
# Activate parallel backend via "parallel" package. (Optional). Uncomment if want to run in sequential mode.
cl <- makeCluster(10)

# Here, we focused on a specific scenario for illustration.
set.seed(3627)
dat <- generateCountData(
  n = 300, p = 1000, K = 2, param = 1, sdsignal = 1.5,
  tag.samples = TRUE, allZero.rm = TRUE, nonzero_prop = .10
)

dds <- createDDSobject(dat)

# Step 2: Filter data (Remove low quality genes with near zero variances.)
dds_processed <- filterCounts(dds)

# 'nFeat' is the number of DE features selected from complete dataset.
nFeat <- 300

# DE Analysis via DESeq2. Select DE features via DESeq2. Number of features is 'nFeat'.
dds_diffExp <- diffExp(dds_processed)
selectedGenes_DESeq <- selectDEfeatures_DESeq(dds_diffExp, nFeatures = nFeat)

# Step 3: Transform filtered counts using VST method available in DESeq2.
# VST transformed counts will be used for ROC-based DE analysis.
dds_processed <- preProcessCounts(dds_processed, normalize = TRUE, transform = TRUE, transformationMethod = "vst")

# Extract sample information from DESeqDataSet object.
col_data <- colData(dds_processed$DESeqObject) |>
  as.data.frame() |>
  mutate(
    condition01 = if_else(condition == "C2", 1, 0),
    condition = factor(condition, levels = c("C1", "C2"))
  ) |>
  DataFrame()

# Prepare data for ROC analysis.
rocData <- bind_cols(
  as_tibble(t(dds_processed$transformedCounts)),
  tibble(response = col_data$condition01)
)

# DE Analysis via ROC-based methods.
# For computational efficiency, we perform parallel processing here.
# If 'cluster' is NULL or not provided, computations will be done sequentially.
dds_diffExp_ROC <- diffExp_ROC(.object = dds_processed, cluster = cl)

# Select DE features via ROC-based methods. Number of features is 'nFeat'.
selectedGenes_ROC <- selectDEfeatures_ROC(dds_diffExp_ROC, nFeat = nFeat)

# Combine selected genes from both methods (i.e., DESeq2 and ROC-based).
selectedGenes <- bind_cols(
  selectedGenes_ROC$selectedFeatures,
  selectedGenes_DESeq$selectedFeatures
)

# Results for each method with the values of metrics.
# Metrics are AUC values for ROC-based methods and log-fold changes for DESeq2.
results <- bind_rows(
  selectedGenes_DESeq$results,
  selectedGenes_ROC$results
)

# Combine all results in a list.
results_DiffExp <- list(
  results = results,
  selectedFeatures = selectedGenes,
  rocData = rocData
)

# True Detection Rates for Improper Genes
improperGenes <- dds_diffExp_ROC$improper_Genes
DEGs <- dds_diffExp_ROC$DE_Genes

TDR_IGs <- tibble(
  LROC = sum(selectedGenes$LROC %in% improperGenes) / length(improperGenes),
  gAUC = sum(selectedGenes$gAUC %in% improperGenes) / length(improperGenes),
  AUC = sum(selectedGenes$AUC %in% improperGenes) / length(improperGenes),
  DESeq = sum(selectedGenes$DESeq %in% improperGenes) / length(improperGenes)
)

DE_Results <- results_DiffExp$results |>
  arrange(Method) |>
  pivot_wider(id_cols = Gene, names_from = Method, values_from = Value) |>
  mutate(
    improper = factor(if_else(Gene %in% improperGenes, 1, 0)),
    DEGs = factor(if_else(Gene %in% DEGs, 1, 0)),
    type = factor(
      case_when(
        (DEGs == 1 & improper == 1) ~ "IG",
        (DEGs == 1 & improper == 0) ~ "DEG",
        .default = "Not significant"
      )
    )
  )

axisTheme <- theme(
  axis.text.x = element_text(margin = margin(t = 5, b = 5)),
  axis.text.y = element_text(margin = margin(r = 5, l = 5))
)
baseFontSize <- 12

LROC_AUC <- ggplot(DE_Results, aes(x = AUC, y = LROC, colour = type)) +
  geom_point(size = 2) +
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_markdown(lineheight = 3)
  ) +
  axisTheme +
  scale_colour_manual(values = c(rgb(1, 0, 0, .2), rgb(0, 0.3, 1, .3), "#0000004D"), labels = c("Diff. Exp. Genes", "Improper Genes", "Not significant")) +
  guides(colour = "none") +
  labs(
    x = "Area Under Curve (cAUC)<br><span style='display:block; text-align:center; font-size:18pt;'>(a)</span>",
    y = "Length of ROC Curve (LROC)"
  )
# guides(colour = guide_legend(override.aes = list(colour = c(rgb(1, 0, 0), rgb(0, 0.3, 1), "#000000")), position = "top", title = NULL))

gAUC_AUC <- ggplot(DE_Results, aes(x = AUC, y = gAUC, colour = type)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, colour = "gray30", linetype = 2) +
  geom_hline(yintercept = .5, linetype = 2, colour = "gray60") +
  geom_vline(xintercept = .5, linetype = 2, colour = "gray60") +
  geom_label(inherit.aes = FALSE, x = .65, y = .45, label = "gAUC < 0.5") +
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(margin = margin(l = 10)),
    axis.title.x = element_markdown(lineheight = 3)
  ) +
  axisTheme +
  scale_colour_manual(values = c(rgb(1, 0, 0, .2), rgb(0, 0.3, 1, .3), "#0000004D"), labels = c("Diff. Exp. Genes", "Improper Genes", "Not significant")) +
  guides(colour = guide_legend(override.aes = list(colour = c(rgb(1, 0, 0), rgb(0, 0.3, 1), "#000000")), position = "top", title = NULL)) +
  labs(
    x = "Area Under Curve (cAUC)<br><span style='display:block; text-align:center; font-size:18pt;'>(b)</span>",
    y = "Generalized Area Under Curve (gAUC)"
  )

gAUC_LROC <- ggplot(DE_Results, aes(x = LROC, y = gAUC, colour = type)) +
  geom_point(size = 2) +
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(margin = margin(l = 10)),
    axis.title.x = element_markdown(lineheight = 3)
  ) +
  axisTheme +
  scale_colour_manual(values = c(rgb(1, 0, 0, .2), rgb(0, 0.3, 1, .3), "#0000004D"), labels = c("Diff. Exp. Genes", "Improper Genes", "Not significant")) +
  guides(colour = "none") +
  labs(
    y = "Generalized Area Under Curve (gAUC)",
    # Başlığı iki satır yapıyoruz; (a) için ayrı stil verebilirsin
    x = "Length of ROC Curve (LROC)<br><span style='display:block; text-align:center; font-size:18pt'>(c)</span>"
  )
# guides(colour = guide_legend(override.aes = list(colour = c(rgb(1, 0, 0), rgb(0, 0.3, 1), "#000000")), position = "top", title = NULL))
gAUC_LROC

smoothData <- DE_Results |>
  filter(type != "Not significant")

gAUC_LROC <- gAUC_LROC +
  geom_smooth(data = smoothData, mapping = aes(x = LROC, y = gAUC, colour = type), se = FALSE) +
  geom_hline(yintercept = .5, linetype = 2, colour = "gray60")

figALL <- (LROC_AUC + gAUC_AUC + gAUC_LROC) +
  plot_layout() &
  theme(plot.margin = margin(t = 0, b = 0, r = 5, l = 5))

# Save the figure
# ggsave(filename = "figure/Simulation/IG_visual_inspect.png", plot = figALL, device = png,
#        width = 28, height = 12, units = "cm", dpi = 320)

# ggsave(filename = "document/manuscript/figure/IG_visual_inspect.png", plot = figALL, device = png,
#        width = 28, height = 12, units = "cm", dpi = 320)
