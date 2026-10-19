REPO_ROOT <- Sys.getenv("BCAT1_REPO", ".")  # publication: set BCAT1_REPO or run from the repository root
#!/usr/bin/env Rscript
# =============================================================================
# PHASE 3 — replicate-blocked differential occupancy, spike-in as an offset.
#
# Design:  ~ replicate + condition
#
# REPLICATE IS A BLOCKING FACTOR, NOT A NUISANCE TO AVERAGE OVER. For H3K27ac the
# two replicates ARE the two batches (QC report 5, 6d: promoter enrichment 3.34x
# vs 6.01x), so blocking is what keeps the batch difference out of the condition
# effect. It costs a residual degree of freedom and will widen the intervals --
# that is the correct response to the data, not something to tune away.
#
# SPIKE-IN AS AN OFFSET, NOT A POST-HOC MULTIPLICATION. lib.size is set to
# genome_depth / spikein, so the modelled abundance is count/depth*spikein while
# the dispersion is still estimated from raw counts. Multiplying the fitted
# values afterwards would rescale the estimate but leave the mean-variance
# relationship describing unscaled data -- same logFC, dishonest uncertainty.
#
# norm.factors are FORCED TO 1. Running calcNormFactors (TMM) here would
# re-estimate composition from the data and silently undo the spike-in, which is
# the classic way spike-in normalisation gets lost.
#
# Usage: Rscript phase3_edger.R <mark> [<mark> ...]
# =============================================================================
suppressPackageStartupMessages({library(edgeR)})

QC  <- file.path(REPO_ROOT, "CUTRUN/analysis_QC_normalization/multimark_2026-07")
DOM <- file.path(QC, "domains")
OUT <- file.path(QC, "data")
FDR_CUT <- 0.05

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- c("H3K4me3", "H3K27ac", "H3K9pr", "ENL")

summary_rows <- list()

for (mark in args) {
  cf <- file.path(DOM, paste0(mark, ".counts.tsv"))
  sf <- file.path(DOM, paste0(mark, ".samples.tsv"))
  if (!file.exists(cf) || !file.exists(sf)) {
    message("[skip] ", mark, ": counts/samples not found"); next
  }
  cnt <- read.delim(cf, check.names = FALSE)
  smp <- read.delim(sf, stringsAsFactors = FALSE)

  ann <- cnt[, c("chrom", "start", "end", "name")]
  mat <- as.matrix(cnt[, smp$label, drop = FALSE])
  rownames(mat) <- ann$name

  # condition: shNT is the reference, so logFC is shBCAT1 relative to shNT
  smp$condition <- factor(smp$condition, levels = c("shNT", "shBCAT1"))
  smp$replicate <- factor(smp$replicate)

  y <- DGEList(counts = mat, samples = smp)
  y$samples$lib.size    <- smp$edger_lib_size     # spike-in enters HERE
  y$samples$norm.factors <- 1                     # never TMM on spike-in data

  keep <- filterByExpr(y, group = smp$condition)
  y <- y[keep, , keep.lib.sizes = TRUE]           # keep the offsets we just set
  ann <- ann[keep, ]

  design <- model.matrix(~ replicate + condition, data = smp)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- glmQLFTest(fit, coef = "conditionshBCAT1")
  tt  <- topTags(res, n = Inf, sort.by = "none")$table

  outdf <- cbind(ann, tt)
  outdf <- outdf[order(outdf$FDR, -abs(outdf$logFC)), ]
  write.csv(outdf, file.path(OUT, paste0(mark, "_domain_differential.csv")),
            row.names = FALSE)

  nsig  <- sum(tt$FDR < FDR_CUT)
  ndown <- sum(tt$FDR < FDR_CUT & tt$logFC < 0)
  nup   <- sum(tt$FDR < FDR_CUT & tt$logFC > 0)
  bcv   <- sqrt(y$common.dispersion)
  cat(sprintf(
    "\n[%s] %d domains tested (of %d)\n  BCV %.3f | residual df %d\n  median logFC %+.3f  (spike-in x%.2f applied as offset)\n  FDR<%.2f: %d  (down %d / up %d)\n",
    mark, nrow(tt), length(keep), bcv, nrow(smp) - ncol(design),
    median(tt$logFC), smp$spikein[smp$condition == "shBCAT1"][1],
    FDR_CUT, nsig, ndown, nup))

  summary_rows[[mark]] <- data.frame(
    mark = mark, domains_tested = nrow(tt), domains_total = length(keep),
    BCV = round(bcv, 4), median_logFC = round(median(tt$logFC), 4),
    n_FDR05 = nsig, n_down = ndown, n_up = nup,
    spikein = smp$spikein[smp$condition == "shBCAT1"][1])
}

if (length(summary_rows)) {
  s <- do.call(rbind, summary_rows)
  write.csv(s, file.path(OUT, "domain_differential_summary.csv"), row.names = FALSE)
  cat("\n=== summary ===\n"); print(s, row.names = FALSE)
  cat("\n[ok] per-mark tables: ", OUT, "/<mark>_domain_differential.csv\n", sep = "")
}
