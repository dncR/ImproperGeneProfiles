library(DESeq2)
library(tidyverse)
library(parallel)
library(ggrepel)
source("R/helper_functions.R")

cl <- makeCluster(10)
saveResults <- FALSE

# 2. Read data from local file ----
count <- read.csv(file.path("data/Cervical", "count.csv"), row.names = 1)
col_data <- read.csv(file.path("data/Cervical", "condition.csv"), row.names = 1)

col_data <- col_data %>%
  mutate(
    condition01 = if_else(condition == "Tumor", 1, 0),
    condition = factor(condition, levels = c("Normal", "Tumor"))
  )

count <- select(count, all_of(col_data$colID))
rownames(col_data) <- col_data$colID

dds <- list(
  DESeqObject = DESeqDataSetFromMatrix(count, col_data, ~condition),
  DE_Genes = NULL,
  improper_Genes = NULL
)

class(dds) <- "dds_raw"

dds_processed <- filterCounts(dds)
nFeat <- nrow(dds_processed$DESeqObject)

# DE Analysis via DESeq2
dds_diffExp <- diffExp(dds_processed, nonzero = TRUE)
selectedGenes_DESeq <- selectDEfeatures_DESeq(dds_diffExp, nFeatures = nFeat)

# Step 3: Transform filtered counts using VST method.
dds_processed <- preProcessCounts(dds_processed,
  normalize = TRUE, transform = TRUE,
  transformationMethod = "vst", nonzero = TRUE
)

rocData <- bind_cols(
  as_tibble(t(dds_processed$transformedCounts)),
  tibble(response = colData(dds_processed$DESeqObject)$condition01)
)

# DE Analysis via ROC-based methods.
dds_diffExp_ROC <- diffExp_ROC(.object = dds_processed, cluster = cl)
selectedGenes_ROC <- selectDEfeatures_ROC(dds_diffExp_ROC, nFeat = nFeat)

selectedGenes <- bind_cols(
  selectedGenes_ROC$selectedFeatures,
  selectedGenes_DESeq$selectedFeatures
)

results <- bind_rows(
  selectedGenes_DESeq$results,
  selectedGenes_ROC$results
)

results_DiffExp <- list(
  results = results,
  selectedFeatures = selectedGenes,
  rocData = rocData,
  DESeq_results = dds_diffExp$results
)

# Save DE Analysis results
if (saveResults) {
  save(results_DiffExp, file = file.path("saved/Cervical", "Cervical_results.Rda"))
}

# 4. Feature Select ----
load(file.path("saved/Cervical", "Cervical_results.Rda"))
geneNames <- rownames(results_DiffExp$DESeq_results)

alpha_threshold <- 0.1
logFC_threshold <- 0.6

DEGs_DESeq <- tibble(Genes = geneNames) |>
  bind_cols(as_tibble(results_DiffExp$DESeq_results)) |>
  mutate(
    log10p = -log10(padj),
    DE_Status01 = case_when(
      abs(log2FoldChange) > logFC_threshold & log10p > -log10(alpha_threshold) ~ 1,
      .default = 0
    )
  ) |>
  filter(DE_Status01 == 1) |>
  pull(Genes)

## 4.1. IGs inspection ----
axisTheme <- theme(
  axis.text.x = element_text(margin = margin(t = 5, b = 5)),
  axis.text.y = element_text(margin = margin(r = 5, l = 5))
)
baseFontSize <- 12

colour_DE <- rgb(1, 0, 0)
colour_IG <- rgb(0, 0.3, 1)
colour_NS <- "#000000"
colour_DE_alpha <- rgb(1, 0, 0, .2)
colour_IG_alpha <- rgb(0, 0.3, 1, .3)
colour_NS_alpha <- "#0000004D"

pointColours <- c(colour_DE_alpha, colour_IG_alpha, colour_NS_alpha)

lroc_threshold <- 1.48
auc_threshold <- .6
gauc_threshold <- .7

DE_Results <- results_DiffExp$results |>
  arrange(Method) |>
  pivot_wider(id_cols = Gene, names_from = Method, values_from = Value) |>
  mutate(
    IGs_LROC = factor(if_else(
      AUC < auc_threshold & LROC > lroc_threshold, 1, 0
    )),
    IGs_gAUC = factor(if_else(
      AUC < auc_threshold & gAUC > gauc_threshold, 1, 0
    ))
  )

### LROC ----
plotData <- DE_Results |>
  mutate(
    type = factor(
      case_when(
        .default = "Not significant",
        IGs_LROC == 1 ~ "Improper Genes",
        Gene %in% DEGs_DESeq & IGs_LROC == 0 ~ "Diff. Exp. Genes"
      )
    )
  )

IGs <- plotData |>
  filter(IGs_LROC == 1) |>
  arrange(desc(LROC)) |>
  pull(Gene)

IGs_unique <- setdiff(IGs, DEGs_DESeq)
if (length(IGs_unique) > 0) {
  IGs_unique <- head(IGs_unique, n = 5)
}

repel_text_data <- plotData |>
  mutate(
    repel_text = case_when(
      .default = NA_character_,
      Gene %in% IGs_unique ~ Gene
    )
  ) |>
  filter(!is.na(repel_text))

LROC_AUC <- ggplot(plotData, aes(x = AUC, y = LROC, colour = type)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = repel_text_data,
    mapping = aes(label = repel_text),
    colour = colour_IG, max.overlaps = 50
  ) +
  geom_point(
    data = repel_text_data,
    color = colour_IG, size = 2
  ) +
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_markdown(lineheight = 3)
  ) +
  axisTheme +
  scale_colour_manual(values = pointColours, labels = c("Diff. Exp. miRNAs", "Improper miRNAs", "Not significant")) +
  guides(colour = "none") +
  labs(
    x = "Area Under Curve (cAUC)<br><span style='display:block; text-align:center; font-size:18pt;'>(a)</span>",
    y = "Length of ROC Curve (LROC)"
  )

###  gAUC ----
plotData <- DE_Results |>
  mutate(
    type = factor(
      case_when(
        .default = "Not significant",
        IGs_gAUC == 1 ~ "Improper miRNAs",
        Gene %in% DEGs_DESeq & IGs_gAUC == 0 ~ "Diff. Exp. miRNAs"
      )
    )
  )

IGs <- plotData |>
  filter(IGs_gAUC == 1) |>
  arrange(desc(gAUC)) |>
  pull(Gene)

IGs_unique <- setdiff(IGs, DEGs_DESeq)
if (length(IGs_unique) > 0) {
  IGs_unique <- head(IGs_unique, n = 5)
}

repel_text_data <- plotData |>
  mutate(
    repel_text = case_when(
      .default = NA_character_,
      Gene %in% IGs_unique ~ Gene
    )
  ) |>
  filter(!is.na(repel_text))


gAUC_AUC <- ggplot(plotData, aes(x = AUC, y = gAUC, colour = type)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = repel_text_data,
    mapping = aes(label = repel_text),
    colour = colour_IG, max.overlaps = 50
  ) +
  geom_point(
    data = repel_text_data,
    color = colour_IG, size = 2
  ) +
  geom_abline(slope = 1, intercept = 0, colour = "gray30", linetype = 2) +
  geom_vline(xintercept = .5, linetype = 2, colour = "gray60") +
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(margin = margin(l = 20)),
    axis.title.x = element_markdown(lineheight = 3)
  ) +
  axisTheme +
  scale_colour_manual(values = c(rgb(1, 0, 0, .2), rgb(0, 0.3, 1, .3), "#0000004D"), labels = c("Diff. Exp. miRNAs", "Improper miRNAs", "Not significant")) +
  guides(colour = guide_legend(override.aes = list(colour = c(rgb(1, 0, 0), rgb(0, 0.3, 1), "#000000")), position = "top", title = NULL)) +
  labs(
    x = "Area Under Curve (cAUC)<br><span style='display:block; text-align:center; font-size:18pt;'>(b)</span>",
    y = "Generalized Area Under Curve (gAUC)"
  ) +
  lims(y = c(.5, 1), x = c(.5, 1))

### gAUC_LROC ----
plotData <- DE_Results |>
  mutate(
    type = factor(
      case_when(
        .default = "Not significant",
        IGs_LROC == 1 | IGs_gAUC == 1 ~ "Improper miRNAs",
        (Gene %in% DEGs_DESeq) & (IGs_gAUC == 0 | IGs_LROC == 0) ~ "Diff. Exp. miRNAs"
      )
    )
  )

IGs_LROC <- plotData |>
  filter(IGs_LROC == 1) |>
  arrange(desc(LROC)) |>
  pull(Gene) |>
  setdiff(DEGs_DESeq) |>
  head(n = 5)

IGs_gAUC <- plotData |>
  filter(IGs_gAUC == 1) |>
  arrange(desc(gAUC)) |>
  pull(Gene) |>
  setdiff(DEGs_DESeq) |>
  head(n = 5)

IGs_unique <- unique(c(IGs_LROC, IGs_gAUC))

repel_text_data <- plotData |>
  mutate(
    repel_text = case_when(
      .default = NA_character_,
      Gene %in% IGs_unique ~ Gene
    )
  ) |>
  filter(!is.na(repel_text))

gAUC_LROC <- ggplot(plotData, aes(x = LROC, y = gAUC, colour = type)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = repel_text_data,
    mapping = aes(label = repel_text),
    colour = colour_IG, max.overlaps = 50
  ) +
  geom_point(
    data = repel_text_data,
    color = colour_IG, size = 2
  ) +
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(margin = margin(l = 20)),
    axis.title.x = element_markdown(lineheight = 3)
  ) +
  axisTheme +
  scale_colour_manual(values = c(rgb(1, 0, 0, .2), rgb(0, 0.3, 1, .3), "#0000004D"), labels = c("Diff. Exp. miRNAs", "Improper miRNAs", "Not significant")) +
  guides(colour = "none") +
  labs(
    y = "Generalized Area Under Curve (gAUC)",
    x = "Length of ROC Curve (LROC)<br><span style='display:block; text-align:center; font-size:18pt'>(c)</span>"
  )

smoothData <- plotData |>
  filter(type != "Not significant")

gAUC_LROC <- gAUC_LROC +
  geom_smooth(data = smoothData, mapping = aes(x = LROC, y = gAUC, colour = type), se = FALSE)

figALL <- (LROC_AUC + gAUC_AUC + gAUC_LROC) +
  plot_layout() &
  theme(plot.margin = margin(t = 0, b = 0, r = 2, l = 2))

print(figALL)

# ggsave(filename = "figure/Cervical/cervical_IG_visual_inspect.jpeg", plot = figALL, device = jpeg,
#        width = 30, height = 12, units = "cm", dpi = 320, create.dir = TRUE)

# ggsave(filename = "document/manuscript/figure/cervical_IG_visual_inspect.jpeg", plot = figALL, device = jpeg,
#        width = 30, height = 12, units = "cm", dpi = 320)

## 4.2 Volcano Plot ----
# Genes selected via DESeq, log2FC between [-0.6, 0.6] and unadjusted p-value below 0.05
geneNames <- rownames(results_DiffExp$DESeq_results)
DESeq_Results <- tibble(Genes = geneNames) |>
  bind_cols(as_tibble(results_DiffExp$DESeq_results))

DE_Results <- results_DiffExp$results |>
  arrange(Method) |>
  pivot_wider(id_cols = Gene, names_from = Method, values_from = Value) |>
  mutate(
    IGs_LROC = factor(if_else(
      AUC < auc_threshold & LROC > lroc_threshold, 1, 0
    )),
    IGs_gAUC = factor(if_else(
      AUC < auc_threshold & gAUC > gauc_threshold, 1, 0
    ))
  )

DESeq_Results <- DESeq_Results |>
  mutate(
    log10p = -log10(pvalue),
    DE_Status = factor(
      case_when(
        log2FoldChange < -logFC_threshold & log10p > -log10(alpha_threshold) ~ -1,
        log2FoldChange > logFC_threshold & log10p > -log10(alpha_threshold) ~ 1,
        .default = 0
      ),
      levels = c(-1, 0, 1),
      labels = c("Downregulated", "Not significant", "Upregulated")
    ),
    DE_Status2 = factor(
      case_when(
        log2FoldChange < -logFC_threshold & pvalue < alpha_threshold ~ -1,
        log2FoldChange > logFC_threshold & pvalue < alpha_threshold ~ 1,
        .default = 0
      ),
      levels = c(-1, 0, 1),
      labels = c("Downregulated", "Not significant", "Upregulated")
    )
  )

DEGs_DESeq <- DESeq_Results |>
  filter(DE_Status != "Not significant") |>
  pull(Genes)

plotData <- DE_Results |>
  mutate(
    type = factor(
      case_when(
        .default = "Not significant",
        IGs_LROC == 1 | IGs_gAUC == 1 ~ "Improper Genes",
        (Gene %in% DEGs_DESeq) & (IGs_gAUC == 0 | IGs_LROC == 0) ~ "Diff. Exp. Genes"
      )
    )
  )

IGs_LROC <- plotData |>
  filter(IGs_LROC == 1) |>
  arrange(desc(LROC)) |>
  pull(Gene) |>
  (\(x){
    remove <- grepl(pattern = "Candidate", x, ignore.case = TRUE)
    x[!remove]
  })()

IGs_gAUC <- plotData |>
  filter(IGs_gAUC == 1) |>
  arrange(desc(gAUC)) |>
  pull(Gene) |>
  (\(x){
    remove <- grepl(pattern = "Candidate", x, ignore.case = TRUE)
    x[!remove]
  })()

IGs <- unique(c(IGs_LROC, IGs_gAUC))
IGs_not_in_DESeq <- setdiff(IGs, DEGs_DESeq)

topIGs_not_in_DESeq <- plotData |>
  filter(Gene %in% IGs_not_in_DESeq) |>
  arrange(desc(gAUC)) |>
  pull(Gene)

IGs_in_DESeq <- IGs[IGs %in% DEGs_DESeq]

# Individual Plots for Improper Genes
topIGs_not_in_DESeq_selected <- head(topIGs_not_in_DESeq, n = 12)
validated_miRNAs <- c(
  "let-7d", "miR-30c-2*", "miR-382", "miR-145", "miR-301a", "miR-19a",
  "miR-19b", "miR-106b", "miR-15a", "miR-598", "miR-32", "miR-454"
)
lapply(topIGs_not_in_DESeq_selected, function(x, .data = results_DiffExp$rocData, fig_path = file.path("figure", "Cervical", "overall"), ...) {
  if (!dir.exists(fig_path)) {
    dir.create(fig_path, recursive = TRUE)
  }

  png(filename = file.path(fig_path, paste0(x, ".png")))
  plot(pROC::roc(predictor = .data[[x]], response = .data$response), main = x)
  dev.off()

  .data$response <- factor(.data$response)
  fig <- ggplot(.data, aes(x = !!sym(x), y = after_stat(density), group = response, colour = response)) +
    geom_density() +
    theme_bw()

  ggsave(filename = file.path(fig_path, paste0(x, "-dens.png")), device = png)
})

jpeg(
  filename = file.path(file.path("figure", "Cervical", "overall"), paste0("cervical_IG_genes_panel_density.jpeg")),
  width = 20, height = 12, units = "cm", res = 450, pointsize = 10
)
par(mfrow = c(3, 4), oma = c(0, 0, 0, 0))
for (i in topIGs_not_in_DESeq_selected) {
  dens1 <- density(results_DiffExp$rocData[[i]][results_DiffExp$rocData$response == 1])
  dens0 <- density(results_DiffExp$rocData[[i]][results_DiffExp$rocData$response != 1])

  xmin <- min(min(dens1$x), min(dens0$x))
  xmax <- max(max(dens1$x), max(dens0$x))
  ymin <- min(min(dens1$y), min(dens0$y))
  ymax <- max(max(dens1$y), max(dens0$y))

  par(mar = c(3, 4, 3, 2))
  plot(dens1, xlim = c(xmin, xmax), ylim = c(ymin, ymax), main = i, xlab = "", ylab = "Density")
  lines(dens0, lty = 2)
}
dev.off()


# Combined plot for Improper Genes (2 rows 4 columns)
jpeg(
  filename = file.path(file.path("figure", "Cervical", "overall"), paste0("cervical_IG_genes_panel.jpeg")),
  width = 20, height = 15, units = "cm", res = 450, pointsize = 10
)
par(mfrow = c(3, 4))
for (i in topIGs_not_in_DESeq_selected) {
  plot(pROC::roc(predictor = results_DiffExp$rocData[[i]], response = results_DiffExp$rocData$response), main = i, lwd = 1)
}
dev.off()


# Biostatsquid theme
theme_set(
  theme_classic(base_size = 12) +
    theme(
      axis.title.y = element_text(face = "bold", margin = margin(0, 10, 0, 0), size = rel(1.1), color = "black"),
      axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(10, 0, 0, 0), size = rel(1.1), color = "black"),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5)),
      plot.margin = unit(c(.1, .1, 0, .1), "cm")
    )
)

plotData <- DESeq_Results

tmp1 <- DESeq_Results |>
  filter(Genes %in% topIGs_not_in_DESeq[1:10])

tmp2 <- DESeq_Results |>
  filter(Genes %in% IGs_in_DESeq) |>
  arrange(desc(abs(log2FoldChange))) |>
  slice_head(n = 6)

repelTextData_volcanoPlot <- bind_rows(tmp1, tmp2)

set.seed(1)
volcanoPlot <- ggplot(data = plotData, aes(x = log2FoldChange, y = log10p, colour = DE_Status2, label = Genes)) +
  geom_vline(xintercept = c(-0.6, 0.6), colour = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), colour = "gray", linetype = "dashed") +
  geom_point(size = 2, position = position_jitter(height = .4)) +
  scale_color_manual(
    values = c(rgb(.3, .8, .3, .2), rgb(.8, .8, .8, .2), rgb(.9, .1, .1, .2)), # to set the colours of our variable
    labels = c("Downregulated", "Not significant", "Upregulated")
  ) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 25)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(x = expression("log"[2] * "FoldChange"), y = expression("-log"[10] * " p-value")) +
  guides(colour = guide_legend(
    override.aes = list(colour = c(rgb(.3, .8, .3), rgb(.8, .8, .8), rgb(.9, .1, .1))),
    position = "top", title = NULL
  )) +
  geom_point(
    mapping = aes(x = log2FoldChange, y = -log10(pvalue)), data = repelTextData_volcanoPlot,
    color = colour_IG, size = 2
  ) +
  geom_text_repel(
    mapping = aes(x = log2FoldChange, y = -log10(pvalue)), data = repelTextData_volcanoPlot,
    color = colour_IG, max.overlaps = Inf
  )

# ggsave(filename = "figure/Cervical/volcanoPlot_cervical.jpeg", plot = volcanoPlot, device = jpeg,
#        width = 18, height = 14, units = "cm", dpi = 320)


# Dispersion estimates for Cervical cancer dataset
count <- read.csv(file.path("data/Cervical", "count.csv"), row.names = 1)
col_data <- read.csv(file.path("data/Cervical", "condition.csv"), row.names = 1)

col_data <- col_data %>%
  mutate(
    condition01 = if_else(condition == "Tumor", 1, 0),
    condition = factor(condition, levels = c("Normal", "Tumor"))
  )

count <- select(count, all_of(col_data$colID))

rownames(col_data) <- col_data$colID
dds <- list(
  DESeqObject = DESeqDataSetFromMatrix(count, col_data, ~condition),
  DE_Genes = NULL,
  improper_Genes = NULL
)
dds_processed <- filterCounts(dds)

dds_disp <- estimateDispersions(estimateSizeFactors(dds_processed$DESeqObject))

disp_plot <- tibble(dispersions = dispersions(dds_disp)) |>
  ggplot(aes(x = dispersions, y = after_stat(density))) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_density(colour = "red") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5, b = 5)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = "Gene-wise dispersion estimates", y = NULL) +
  lims(x = c(0, 10))

ggsave(filename = "cervical_dispersion_estimates.jpeg", plot = disp_plot, device = jpeg, width = 18, height = 14, units = "cm", dpi = 320, path = "figure/Cervical")
