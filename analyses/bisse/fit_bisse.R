#!/usr/bin/env Rscript
# BiSSE fits over the case-study replicates.
#
# For every replicate tree, and for each of the two characters (the driver family's
# presence and the control family's presence), fit the full six-parameter BiSSE model
# and the state-independent constraint (lambda1 ~ lambda0, mu1 ~ mu0), and report the
# likelihood-ratio test: chisq = 2 * (lnL_full - lnL_null), df = 2.
#
# This is deliberately the textbook recipe (diversitree's own vignette): make.bisse,
# starting.point.bisse, find.mle, constrain — defaults throughout (root=ROOT.OBS,
# condition.surv=TRUE, sampling.f=1), because the question is what the method as
# people actually run it reports on these data. Invariant characters (all tips 0 or
# all tips 1) are skipped, as they would be in practice.
#
#   conda run -n bisse Rscript fit_bisse.R [n_cores]

suppressPackageStartupMessages({
  library(diversitree)
  library(ape)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[[1]]) else 8L

file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
HERE <- normalizePath(dirname(sub("--file=", "", file_arg)))
tree_files <- sort(list.files(file.path(HERE, "data", "trees"),
                              pattern = "\\.nwk$", full.names = TRUE))
cat("trees:", length(tree_files), " cores:", n_cores, "\n")

fit_one <- function(nwk_path) {
  tag <- sub("\\.nwk$", "", basename(nwk_path))
  tr <- read.tree(nwk_path)
  st <- read.delim(file.path(HERE, "data", "states", paste0(tag, ".tsv")))
  out <- list()
  for (char in c("driver", "control")) {
    x <- setNames(as.integer(st[[char]]), st$tip)[tr$tip.label]
    k <- sum(x); n <- length(x)
    row <- data.frame(tag = tag, character = char, n_tips = n, n_present = k,
                      lnL_full = NA_real_, lnL_null = NA_real_,
                      lambda0 = NA_real_, lambda1 = NA_real_,
                      mu0 = NA_real_, mu1 = NA_real_,
                      q01 = NA_real_, q10 = NA_real_,
                      chisq = NA_real_, p = NA_real_,
                      status = "skipped_invariant", stringsAsFactors = FALSE)
    if (k > 0 && k < n) {
      res <- try({
        lik  <- make.bisse(tr, x)
        p0   <- starting.point.bisse(tr)
        fit  <- find.mle(lik, p0)
        lik0 <- constrain(lik, lambda1 ~ lambda0, mu1 ~ mu0)
        fit0 <- find.mle(lik0, p0[argnames(lik0)])
        chi  <- 2 * (fit$lnLik - fit0$lnLik)
        row$lnL_full <- fit$lnLik
        row$lnL_null <- fit0$lnLik
        cf <- coef(fit)
        row$lambda0 <- cf[["lambda0"]]; row$lambda1 <- cf[["lambda1"]]
        row$mu0 <- cf[["mu0"]]; row$mu1 <- cf[["mu1"]]
        row$q01 <- cf[["q01"]]; row$q10 <- cf[["q10"]]
        row$chisq <- chi
        row$p <- pchisq(max(chi, 0), df = 2, lower.tail = FALSE)
        row$status <- "ok"
      }, silent = TRUE)
      if (inherits(res, "try-error"))
        row$status <- paste0("error: ", gsub("[\t\n]", " ", substr(as.character(res), 1, 120)))
    }
    out[[char]] <- row
  }
  do.call(rbind, out)
}

t0 <- Sys.time()
rows <- mclapply(tree_files, fit_one, mc.cores = n_cores)
bad <- sum(vapply(rows, function(r) inherits(r, "try-error") || is.null(r), logical(1)))
if (bad > 0) cat("WARNING:", bad, "replicates returned nothing\n")
res <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
write.table(res, file.path(HERE, "fits.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("wrote", file.path(HERE, "fits.tsv"), "-", nrow(res), "rows in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("status counts:\n")
print(table(res$status == "ok"))
