#!/usr/bin/env Rscript
# BiSSE fits for the tree-size arms. Same textbook recipe as fit_bisse.R, but over
# data_size/ and restricted to one size prefix per invocation, so the 500-tip fits
# land before the 1,000-tip fits begin.
#
#   conda run -n bisse Rscript fit_bisse_size.R n500 fits_n500.tsv 8

suppressPackageStartupMessages({
  library(diversitree)
  library(ape)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[[1]]
outfile <- args[[2]]
n_cores <- if (length(args) >= 3) as.integer(args[[3]]) else 8L

file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
HERE <- normalizePath(dirname(sub("--file=", "", file_arg)))
tree_files <- sort(list.files(file.path(HERE, "data_size", "trees"),
                              pattern = paste0("^", prefix, "_.*\\.nwk$"),
                              full.names = TRUE))
cat("prefix:", prefix, " trees:", length(tree_files), " cores:", n_cores, "\n")

fit_one <- function(nwk_path) {
  tag <- sub("\\.nwk$", "", basename(nwk_path))
  tr <- read.tree(nwk_path)
  st <- read.delim(file.path(HERE, "data_size", "states", paste0(tag, ".tsv")))
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
res <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
write.table(res, file.path(HERE, outfile),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("wrote", outfile, "-", nrow(res), "rows in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
print(table(res$status == "ok"))
