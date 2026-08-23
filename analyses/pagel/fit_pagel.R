# Pagel's correlated-evolution test on every replicate: habitat against family A,
# and habitat against the independent control, with phytools::fitPagel.
#
#   Rscript fit_pagel.R          # writes fits.tsv
suppressPackageStartupMessages({library(phytools); library(ape); library(parallel)})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
arms <- c("feedback", "trait2gen", "gen2trait", "null")
jobs <- expand.grid(arm = arms, rep = 1:150, stringsAsFactors = FALSE)

one <- function(i) {
  arm <- jobs$arm[i]; rep <- jobs$rep[i]
  tree <- read.tree(file.path(here, "data", "trees", sprintf("%s_r%03d.nwk", arm, rep)))
  tips <- read.delim(file.path(here, "data", "states", sprintf("%s_r%03d.tsv", arm, rep)))
  rownames(tips) <- tips$tip
  tips <- tips[tree$tip.label, ]
  out <- list()
  for (pair in c("A", "ctrl")) {
    x <- setNames(tips$habitat, tips$tip)
    y <- setNames(tips[[pair]], tips$tip)
    if (length(unique(x)) < 2 || length(unique(y)) < 2) {
      out[[pair]] <- data.frame(arm = arm, rep = rep, pair = pair, P = NA, note = "invariant")
      next
    }
    fit <- tryCatch(fitPagel(tree, x, y), error = function(e) NULL)
    out[[pair]] <- if (is.null(fit))
      data.frame(arm = arm, rep = rep, pair = pair, P = NA, note = "fit-error")
    else
      data.frame(arm = arm, rep = rep, pair = pair, P = fit$P, note = "")
  }
  do.call(rbind, out)
}

res <- mclapply(seq_len(nrow(jobs)), one, mc.cores = 8)
out <- do.call(rbind, res)
write.table(out, file.path(here, "fits.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("fits written:", nrow(out), "\n")
