# Custom functions ----

printStatus <- function(idx, current, N = NULL) {
  cat(rep("\n", 50))
  cat("Running simulation...", "\n", sep = "")
  cat("Progress: ", idx, " out of ", N, " (", round(100 * (idx - 1) / N, 1), "% completed.)", "\n\n", sep = "")
  cat("Simulation Parameters:", "\n")
  cat("  - Simulation ID:", current$simID, "\n")
  cat("  - Sample size (n):", current$n, "\n")
  cat("  - Number of genes (p): fixed at", current$p, "\n")
  cat("  - Prop. of differentially expressed genes (DE):", current$propDE, "\n")
  cat("  - Prop. of improper genes (IR) within DE genes:", current$propIR, "\n")
  cat(
    "  - Overdispersion (phi):", current$phi,
    case_when(
      current$phi == 0.01 ~ "(very slight)",
      current$phi == 0.1 ~ "(moderate)",
      current$phi == 1 ~ "(very high)",
    ), "\n"
  )
  cat("  - Offset parameter (controlling log2-FoldChange):", current$sdsignal, "\n\n")
}

#' @title Generate Count Data from RNA-Sequencing
#'
#' @description
#' Generate raw gene expression data from an RNA-Sequencing experiment. Generated data are mapped read counts, which are related to the expression levels of genes (or gene regions).
#'
#' @param n integer. Number of samples (i.e., columns) in the generated data.
#' @param p integer. Number of features (i.e., rows) in the generated data.
#' @param K integer. Number of classes. Default is 2.
#' @param param numeric. Dispersion parameter. Overdispersion parameter of Negative Binomial distribution is equal to 1/param. Mapped read counts are assumed to be distributed as Negative Binomial with overdispersion parameter \code{1/param}.
#' @param sdsignal The extent to which the classes are different. If this equals zero then there are no class differences and if this is large then the classes are very different.
#' @param DE numeric. Probability of a gene being differentially expressed between classes. It should be a numeric value between 0 and 1. Default is 0.3.
#' @param allZero.rm a logical. If TRUE, features whose values are all zero will be removed from generated dataset. This is not common, yet, it might be observed when sample size 'n' is small. Default is TRUE.
#' @param tag.samples a logical. If TRUE, samples are tagged into one of generated classes, e.g., control vs. healthy.
#'
#' @return a list with all elements of generated data set.
#' @export
#'
#' @examples
#' 1L
generateCountData <- function(n = 20, p = 20, K = 2, param = 10, sdsignal = 1, DE = 0.3,
                              IR = 0.25, nonzero_prop = 0, min_count = 1,
                              allZero.rm = TRUE, tag.samples = FALSE, ...) {
  # n = 500; p = 3000; K = 2; param = 2; sdsignal = 2; DE = 0.3;
  # IR = 0.25; nonzero_prop = .10; allZero.rm = TRUE; tag.samples = FALSE
  # min_count = 1

  add_mincount <- function(.x, .nz_p = .10, .val = 1) {
    zeros_l <- (apply(.x, 1, min) == 0)
    p_nz <- sum(!zeros_l) / length(zeros_l)

    # If the proportion of genes having zero cells is below 10%, then return count
    # matrix withoud adding constant value
    if (p_nz > .10) {
      return(.x)
    }

    index <- as.numeric(
      sample(
        which(zeros_l),
        ceiling(length(zeros_l) * (.nz_p - p_nz)),
        FALSE
      )
    )

    .x[index, ] <- .x[index, ] + .val
    return(.x)
  }

  if (n < 4 * K) {
    stop("We require n to be at least 4*K.")
  }

  if (DE < 0 || DE > 1) {
    stop("DE should be a value within [0, 1].")
  }

  if (IR < 0 || IR > 1) {
    stop("IR should be a value within [0, 1].")
  }

  q0 <- rexp(p, rate = 1 / 25) # Gene total, g.j
  isDE <- runif(p) <= DE
  isImproper <- (isDE & runif(p) <= IR)
  classk <- fc <- matrix(NA, nrow = K + 1, ncol = p)

  for (k in 1:K) {
    # dkj parameter
    lfc <- rnorm(p, mean = 0, sd = sdsignal)

    ## Hangi genlerin DE oldugu bilgisi. BunlarD1 da generate edilen veri ile dC6ndC<rebiliriz.
    # classk[k, ] <- ifelse(isDE, q0 * exp(lfc), q0)
    classk[k, ] <- (1 - as.numeric(isDE)) * q0 + as.numeric(isDE) * q0 * exp(lfc)
  }

  classk[3, ] <- classk[1, ]
  ratio <- classk[2, ] / classk[1, ]
  classk[3, isImproper] <- classk[2, isImproper] * ratio[isImproper]

  truesf <- runif(n) * 2 + 0.2
  truesfte <- runif(n) * 2 + 0.2

  # Two classes are set.
  class_set <- rep(1:K, each = floor(n / K))
  if (n %% 2 == 1) {
    class_set <- c(class_set, sample(class_set, 1))
  }
  conds <- sample(class_set, n, FALSE)
  condste <- sample(class_set, n, FALSE)
  x <- xte <- matrix(NA, nrow = n, ncol = p)

  for (i in 1:n) {
    # If given feature is improper, positive class is defined on left and right side.
    # Controls are in the middle.
    if (conds[i] == 2) {
      idx <- ifelse(runif(1) <= 0.50, 1, 3)
      x[i, ] <- rnbinom(p, mu = truesf[i] * classk[idx, ], size = param)
      xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[idx, ], size = param)
    } else {
      idx <- 2
      x[i, ] <- rnbinom(p, mu = truesf[i] * classk[idx, ], size = param)
      xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[idx, ], size = param)
    }
  }

  # Removes all zero columns, if exists.
  if (allZero.rm) {
    rm <- apply(x, 2, sum) == 0
  } else {
    rm <- logical(ncol(x))
  }

  colnames(x) <- colnames(xte) <- paste("G", 1:p, sep = "")

  if (tag.samples) {
    rownames(x) <- rownames(xte) <- paste("S", 1:n, sep = "")
  }

  count <- t(x[, !rm])
  count_te <- t(xte[, !rm])

  if (nonzero_prop > 0) {
    count <- add_mincount(count, .nz_p = nonzero_prop, .val = min_count)
    count_te <- add_mincount(count_te, .nz_p = nonzero_prop, .val = min_count)
  }

  return(
    structure(
      list(
        x = count, xte = count_te, y = conds, yte = condste,
        truesf = truesf, truesfte = truesfte,
        DE_Genes = rownames(count)[isDE[!rm]],
        improper_Genes = rownames(count)[isImproper[!rm]]
      ),
      class = "count.data"
    )
  )
}


# .object: a list returned from "generateCountData" function.
createDDSobject <- function(.object, ...) {
  cnts <- .object$x
  nClasses <- length(unique(.object$y))
  classes <- paste0("C", 1:nClasses)[.object$y]
  genes <- rownames(cnts)

  coldata <- data.frame(
    condition = factor(classes)
  )

  rownames(coldata) <- colnames(cnts)

  dds <- DESeqDataSetFromMatrix(
    countData = cnts,
    colData = coldata,
    design = ~condition
  )

  return(
    structure(
      list(
        DESeqObject = dds,
        DE_Genes = .object$DE_Genes,
        improper_Genes = .object$improper_Genes
      ),
      class = "dds_raw"
    )
  )
}

# .object: an object returned from "createDDSobject" function.
filterCounts <- function(.object, filterLowCounts = TRUE, filterNearZeroVariance = TRUE,
                         lc.threshold = 10, smallestGroupSize = 3, ...) {
  # .object <- dds2
  # filterLowCounts = TRUE
  # filterNearZeroVariance = TRUE
  # lc.threshold = 10
  # smallestGroupSize = 3

  dds <- .object$DESeqObject

  if (!inherits(dds, "DESeqDataSet")) {
    stop("'object' must be of 'DESeqDataSet' class generated by the function 'DESeqDataSetFromMatrix'.")
  }

  filterRes <- list()

  # Near-zero variance filtering
  if (filterNearZeroVariance) {
    filtered_Genes_nzVar <- NULL
    nzVar <- caret::nearZeroVar(t(counts(dds)))
    if (length(nzVar) > 0) {
      filtered_Genes_nzVar <- rownames(dds)[nzVar]
      dds <- dds[-nzVar, ]
      .object$DE_Genes <- setdiff(.object$DE_Genes, filtered_Genes_nzVar)
      .object$improper_Genes <- setdiff(.object$improper_Genes, filtered_Genes_nzVar)
    }

    filterRes[["nearZeroVar"]] <- filtered_Genes_nzVar
  }

  # Keep features having mapped read counts above "lc.threshold" from at least "smallestGroupSize" samples.
  # Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal
  # number of samples. The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation
  # for the minimal number of samples is to specify the smallest group size, e.g. here there
  # are 3 treated samples.
  if (filterLowCounts) {
    filtered_Genes_lowCounts <- NULL
    keep <- rowSums(counts(dds) >= lc.threshold) >= smallestGroupSize
    if (any(!keep)) {
      filtered_Genes_lowCounts <- rownames(dds)[!keep]
      dds <- dds[keep, ]
      .object$DE_Genes <- setdiff(.object$DE_Genes, filtered_Genes_lowCounts)
      .object$improper_Genes <- setdiff(.object$improper_Genes, filtered_Genes_lowCounts)
    }

    filterRes[["lowCounts"]] <- filtered_Genes_lowCounts
  }

  .object$DESeqObject <- dds
  .object$filterRes <- filterRes
  attr(.object, "class") <- "dds_filtered"

  return(.object)
}


# object: an object returned from "createDDSobject" function or filtered object returned by
# "filterCounts" function.
preProcessCounts <- function(.object, normalize = TRUE, transform = TRUE,
                             transformationMethod = c("vst", "rlog", "logCPM"), nonzero = FALSE, ...) {
  # .object = dds_processed
  # .object = dds_diffExp
  # normalize = TRUE
  # transform = TRUE
  # transformationMethod = "vst"

  dds <- .object$DESeqObject
  if (!(class(.object) %in% c("dds_filtered", "dds_diffExp")) | !inherits(dds, "DESeqDataSet")) {
    stop("'object' must be of 'DESeqDataSet' class generated by the function 'DESeqDataSetFromMatrix'.")
  }


  if (nonzero & min(counts(dds)) == 0) {
    counts(dds) <- counts(dds) + 1L
    warning("Adding constant '1' to raw counts, re-calculating...")
  }

  transformationMethod <- match.arg(transformationMethod)

  normalizedCounts <- vstCounts <- NULL
  if (normalize) {
    success <- try({
      dds <- estimateSizeFactors(dds)
    })

    if (!inherits(success, "try-error")) {
      normalizedCounts <- counts(dds, normalize = TRUE)
    } else {
      normalizedCounts <- NULL
    }
  }

  if (transform) {
    if (transformationMethod == "vst") {
      success <- try({
        dds <- varianceStabilizingTransformation(dds)
      })

      if (!inherits(success, "try-error")) {
        vstCounts <- assay(dds)
      } else {
        vstCounts <- NULL
      }
    }
  }

  # .object$DESeqObject <- dds
  .object$normalizedCounts <- normalizedCounts
  .object$transformedCounts <- vstCounts
  attr(.object, "class") <- "dds_transformed"

  return(.object)
}

diffExp <- function(.object, nonzero = FALSE, ...) {
  dds <- .object$DESeqObject

  if (!inherits(dds, "DESeqDataSet")) {
    stop("'object' must be of 'DESeqDataSet' class generated by the function 'DESeqDataSetFromMatrix'.")
  }

  if (nonzero & min(counts(dds)) == 0) {
    counts(dds) <- counts(dds) + 1L
    warning("Adding constant '1' to raw counts, re-calculating...")
  }

  dds <- DESeq(dds)
  res <- results(dds)

  .object$results <- res
  attr(.object, "class") <- "dds_diffExp"

  return(.object)
}


selectDEfeatures_DESeq <- function(.object, nFeatures = 5, ...) {
  tmp <- .object$results
  tmp <- tmp[order(abs(tmp$padj), decreasing = FALSE), ]
  selectedFeatures <- rownames(tmp)[1:nFeatures]
  rm(tmp)

  res <- tibble(
    Gene = rownames(.object$results),
    Method = "DESeq",
    Value = .object$results$log2FoldChange,
    DE = 0
  ) |>
    mutate(
      DE = if_else(Gene %in% selectedFeatures, 1, 0)
    )

  return(
    structure(
      list(
        results = res,
        selectedFeatures = tibble(DESeq = selectedFeatures),
        DE_Genes = .object$DE_Genes,
        improper_Genes = .object$improper_Genes
      ),
      class = "selectedFeatures_DESeq"
    )
  )
}

# controls: negative class
# cases: positive class
# direction: define the location of the distribution of cases as compared to controls.
# direction_stat: define the central tendency measure to define the center of distributions.
#                 Since LROC is based on parametric binormal ROC curve, "mean" is selected as default
#                 statistic to find the correct direction.
nLhat <- function(controls, cases, direction = c("auto", "higher", "lower"),
                  direction_stat = c("mean", "median"), ...) {
  direction <- match.arg(direction)
  direction_stat <- match.arg(direction_stat)

  stat_fun <- switch(direction_stat,
    mean = mean,
    median = median
  )

  x <- controls
  y <- cases

  if (direction != "higher") {
    stat_controls <- do.call(stat_fun, list(x = controls, na.rm = TRUE))
    stat_cases <- do.call(stat_fun, list(x = cases, na.rm = TRUE))

    if (direction == "lower" | (direction == "auto" & (stat_controls > stat_cases))) {
      x <- cases
      y <- controls
    }
  }

  p <- seq(10e-6, 0.99999, length.out = 10000)
  # CHECK THIS ----
  # Burada median ile bir yön belirlemek istediğimizde formül parametric binormal ROC ile hesaplandığı için
  # LROC değeri ile median'ın belirlediği yön uyumsuz oluyor. delta ve ro hatalı hesaplanıyor. Bunun yerine
  # median üzerinden giden ve aşağıdaki tekniğin non-parametrik varyasyonunu kullanan bir yaklaşım üzerinde
  # düşünülebilir.
  #
  # sigmax = 8
  # sigmay = 8
  # mux = 70
  # xup <- mux + 3 * sigmax
  # muy = xup + sigmay * 3
  # ydown <- muy - 3 * sigmay

  mux <- mean(x)
  muy <- mean(y)
  sigmax <- sd(x)
  sigmay <- sd(y)

  x_threshold <- mux + 3 * sigmax
  y_threshold <- muy - 3 * sigmay

  if (x_threshold <= y_threshold) {
    return(1.99998)
  }

  ro <- sigmax / sigmay
  delta <- (muy - mux) / sigmay
  ROC <- 1 - pnorm(qnorm(1 - p, mux, sigmax), mux + delta * sigmax / ro, sigmax / ro)
  ROCprima <- (ro * exp(-0.5 * (delta + ro * qnorm(p))^2)) / exp(-0.5 * (qnorm(p))^2)
  Laux <- numeric()
  Laux[1] <- ifelse(
    ROCprima[1] <= 1,
    sqrt(1 + ROCprima[1]^2) * p[1],
    sqrt(1 + ROCprima[1]^2) / ROCprima[1] * ROC[1]
  )

  for (k in 2:length(p)) {
    Laux[k] <- ifelse(
      ROCprima[k] <= 1,
      sqrt(1 + ROCprima[k]^2) * (p[k] - p[k - 1]),
      sqrt(1 + ROCprima[k]^2) / (ROCprima[k]) * (ROC[k] - ROC[k - 1])
    )
  }
  Lhat <- sum(Laux)
  Lhat <- pmax(Lhat, sqrt(2))
  Lhat <- ifelse(Lhat > 2, 2, Lhat)
  return(Lhat)
}

diffExp_ROC <- function(.object, cluster = NULL, ...) {
  # packages needed on master; workers will load them explicitly
  library(dplyr)

  if (class(.object) != "dds_transformed") {
    stop("'.object' must be of class 'dds_transformed'.")
  }

  col_data <- colData(.object$DESeqObject) |>
    as_tibble() |>
    mutate(
      condition01 = as.numeric(condition) - 1
    )

  # Prepare data frame where columns are markers (genes) and an added response column
  roc_data <- as_tibble(t(.object$transformedCounts)) |>
    mutate(response = col_data$condition01)

  genes <- rownames(.object$DESeqObject)

  # Helper that computes ROC metrics for a single gene
  .roc_worker_fun <- function(g, .improper_genes = NULL) {
    # TODO: ----
    #  "direction" belirlemek için bir algoritma yazılabilir. Burada her senaryoyu deneyip max olan AUC
    # değerini alıyoruz ancak bu durum computational olarak effective değil. Bunun yerine improper olup olmadığını
    # objective bir yaklaşım ile belirleyen algoritma yazılabilir.
    # Initial function used to find direction of ROC curve from marker and response information
    # find_direction <- function(.marker, .response, .positive = NULL){
    #   observed_levels <- unique(.response)
    #
    #   if (!is.factor(.response)){
    #     .response <- factor(.response)
    #   }
    #
    #   .response_levels <- levels(.response)
    #   num_levels <- length(.response_levels)
    #
    #   if (is.null(.positive)){
    #     if (any(c(num_levels, length(observed_levels)) != 2)){
    #       stop("Exactly two levels must be observed for 'response' variable.")
    #     }
    #     .positive <- .response_levels[2]
    #   } else {
    #     if (as.character(.positive) %in% .response_levels){
    #       stop("'positive' class for response variable is not a valid category.")
    #     }
    #   }
    #
    #   .marker_split <- split(.marker, .response)
    #
    #   idx_pos <- which(names(.marker_split) == .positive)
    #   idx_neg <- which(names(.marker_split) != .positive)
    #
    #   stats <- lapply(.marker_split, median) |>
    #     unlist()
    #
    #   ifelse(stats[[idx_pos]] > stats[[idx_neg]], "", "")
    #
    #
    # }

    response <- roc_data[["response"]]
    marker <- roc_data[[g]]

    # LROC
    .data_split <- split(marker, response)
    cases <- .data_split[["1"]]
    controls <- .data_split[["0"]]
    len_roc <- nLhat(controls, cases, direction_stat = "mean", direction = "auto")

    # gAUC
    if (!is.null(.improper_genes)) {
      direction <- if_else(g %in% improper_Genes, "both", "auto")
      gAUC <- nsROC::gROC(marker, response, side = direction)[["auc"]]
    } else {
      gAUC <- max(
        nsROC::gROC(marker, response, side = "auto")[["auc"]],
        nsROC::gROC(marker, response, side = "both")[["auc"]],
        nsROC::gROC(marker, response, side = "both2")[["auc"]]
      )
    }

    # AUC
    AUC <- as.numeric(pROC::auc(response, marker, direction = "auto", quiet = TRUE))
    AUC <- ifelse(AUC < .5, 1 - AUC, AUC)

    tibble(Gene = g) %>%
      dplyr::cross_join(
        tibble(
          Method = c("LROC", "AUC", "gAUC"),
          Value  = c(len_roc, AUC, gAUC)
        )
      )
  }

  # Run either in parallel or serial
  # If 'cluster' is provided and not NULL, use parallel;
  # otherwise (cluster is NULL or missing), run sequentially.
  if (!missing(cluster) && !is.null(cluster)) {
    # Make 'improper_Genes' visible on workers
    improper_Genes <- .object$improper_Genes

    # Export needed objects and (optionally) user-defined functions to workers
    parallel::clusterExport(
      cl = cluster,
      varlist = c("roc_data", "improper_Genes", ".roc_worker_fun", "nLhat"),
      envir = environment()
    )

    # Ensure required packages are loaded on workers
    parallel::clusterEvalQ(cluster, {
      library(dplyr)
      library(tibble)
      library(pROC)
      library(nsROC)
      NULL
    })

    rocRes <- parallel::parLapply(cl = cluster, X = genes, fun = .roc_worker_fun, .improper_genes = improper_Genes)
  } else {
    improper_Genes <- .object$improper_Genes
    rocRes <- lapply(genes, .roc_worker_fun)
  }

  rocRes <- dplyr::bind_rows(rocRes)

  structure(
    list(
      DESeqObject    = .object$DESeqObject,
      resultsDE_ROC  = rocRes,
      DE_Genes       = .object$DE_Genes,
      improper_Genes = .object$improper_Genes,
      nGenes         = nrow(.object$DESeqObject)
    ),
    class = "diffExp_ROC"
  )
}

selectDEfeatures_ROC <- function(.object, nFeat = 10, ...) {
  if (class(.object) != "diffExp_ROC") {
    stop("'object' must be of class 'diffExp_ROC'.")
  }

  if (nFeat > .object$nGenes) {
    warning("Number of features to be selected (nFeat = ", nFeat, ") is larger than number of available features (i.e., ", object$nGenes, "). \n Returning all features.")
  }
  nFeat <- pmin(nFeat, .object$nGenes)

  feature_list <- .object$resultsDE_ROC %>%
    group_by(Method) %>%
    slice_max(n = nFeat, order_by = Value) |>
    slice_head(n = nFeat)

  selectedFeatures <- tibble(
    LROC = pull(filter(feature_list, Method == "LROC"), "Gene"),
    AUC = pull(filter(feature_list, Method == "AUC"), "Gene"),
    gAUC = pull(filter(feature_list, Method == "gAUC"), "Gene")
  )

  res <- .object$resultsDE_ROC |>
    mutate(
      DE = 0
    ) |>
    mutate(
      DE = if_else(Method == "LROC" & Gene %in% selectedFeatures$LROC, 1, DE),
      DE = if_else(Method == "AUC" & Gene %in% selectedFeatures$AUC, 1, DE),
      DE = if_else(Method == "gAUC" & Gene %in% selectedFeatures$gAUC, 1, DE)
    )

  return(
    structure(
      list(
        results = res,
        selectedFeatures = selectedFeatures,
        DE_Genes = .object$DE_Genes,
        improper_Genes = .object$improper_Genes
      ),
      class = "selectedFeatures_ROC"
    )
  )
}

calculatePerformanceMetrics <- function(objectDESeq, objectROC, .n, ...) {
  # objectDESeq: an object of differential expression analysis results returned from DESeq workflow.
  # objectROC: an object of differential expression analysis results returned from ROC workflow.
  # .n: number of features to be selected.
  # ...: additional arguments to be passed to the function. Currently not used.

  # objectDESeq <- dds_diffExp
  # objectROC <- dds_diffExp_ROC
  # .n <- 150

  selectedGenes_DESeq <- selectDEfeatures_DESeq(objectDESeq, nFeatures = .n)
  selectedGenes_ROC <- selectDEfeatures_ROC(objectROC, nFeat = .n)

  results_DESeq <- selectedGenes_DESeq$results
  results_ROC <- selectedGenes_ROC$results

  results <- bind_rows(
    results_DESeq,
    results_ROC
  ) |>
    mutate(
      DE_actual = 0,
      IG_actual = 0
    )

  idx <- results$Gene %in% objectDESeq$DE_Genes
  results[idx, "DE_actual"] <- 1

  idx <- results$Gene %in% objectDESeq$improper_Genes
  results[idx, "IG_actual"] <- 1

  DE_methods <- c("DESeq", "LROC", "AUC", "gAUC")
  DEGs <- lapply(DE_methods, function(m) {
    results |>
      filter(Method == m) |>
      summarise(
        nGenes_selected = .n,
        Method = m,
        group = "DEGs",
        TPR = sum(DE == 1 & DE_actual == 1) / sum(DE_actual == 1),
        TNR = sum(DE == 0 & DE_actual == 0) / sum(DE_actual == 0),
        PPV = sum(DE == 1 & DE_actual == 1) / sum(DE == 1),
        NPV = sum(DE == 0 & DE_actual == 0) / sum(DE == 0)
      )
  }) |>
    bind_rows()

  IGs <- lapply(DE_methods, function(m) {
    results |>
      filter(Method == m) |>
      summarise(
        nGenes_selected = .n,
        Method = m,
        group = "IGs",
        TPR = sum(DE == 1 & IG_actual == 1) / sum(IG_actual == 1),
        TNR = sum(DE == 0 & IG_actual == 0) / sum(IG_actual == 0),
        PPV = sum(DE == 1 & IG_actual == 1) / sum(DE == 1),
        NPV = sum(DE == 0 & IG_actual == 0) / sum(DE == 0)
      )
  }) |>
    bind_rows()

  # results |>
  #   filter(Method == "DESeq") |>
  #   (\(x) table(x$DE, x$DE_actual))()

  res <- bind_rows(DEGs, IGs)
  return(res)
}
