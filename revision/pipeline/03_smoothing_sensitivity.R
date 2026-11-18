#!/usr/bin/env Rscript
# Smoothing sensitivity of the profile-minima analysis (referee 1, P5 l.35).
#
# Re-uses the row-averaged profiles stored by analyse_profiles.R
# (analysis/profiles.rds) and repeats the minima detection for a range of
# smoothing-spline settings, including no smoothing at all.  For every
# image and setting the mean spacing (in pixels), the number of spacings,
# the equivalent degrees of freedom of the spline and the extreme spacing
# ratios relative to the image MEDIAN spacing are written, so that the
# revised two-sided acceptance rule can be applied downstream.
suppressPackageStartupMessages(library(pracma))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE)
revision_dir <- dirname(dirname(normalizePath(script_path, mustWork = TRUE)))
output_dir <- file.path(revision_dir, "analysis")

profiles <- readRDS(file.path(output_dir, "profiles.rds"))$raw
records <- read.csv(file.path(output_dir, "profiles_per_image.csv"), stringsAsFactors = FALSE)
stopifnot(length(profiles) == nrow(records))

settings <- data.frame(
    name = c("none", "spar_0.01", "spar_0.05", "spar_0.10", "spar_0.20",
             "spar_0.40", "spar_0.70", "GCV"),
    spar = c(NA, 0.01, 0.05, 0.1, 0.2, 0.4, 0.7, NA),
    stringsAsFactors = FALSE
)

measure <- function(profile, name, spar) {
    df <- NA_real_
    if (name == "none") {
        y <- profile
    } else if (name == "GCV") {
        fit <- smooth.spline(seq_along(profile), profile, cv = FALSE)
        y <- fit$y; df <- fit$df
    } else {
        fit <- smooth.spline(seq_along(profile), profile, spar = spar)
        y <- fit$y; df <- fit$df
    }
    peaks <- findpeaks(-y, nups = 1, ndowns = 1, threshold = 0.02)
    positions <- if (is.null(peaks)) numeric() else sort(peaks[, 2])
    spacings <- diff(positions)
    if (length(spacings) < 1) {
        return(data.frame(n_spacings = 0, period_px = NA, min_ratio = NA, max_ratio = NA,
                          min_over_mean = NA, max_over_mean = NA, df = df))
    }
    med <- median(spacings)
    data.frame(n_spacings = length(spacings), period_px = mean(spacings),
               min_ratio = min(spacings) / med, max_ratio = max(spacings) / med,
               min_over_mean = min(spacings) / mean(spacings),
               max_over_mean = max(spacings) / mean(spacings), df = df)
}

out <- vector("list", length(profiles) * nrow(settings))
k <- 0L
for (i in seq_along(profiles)) {
    for (j in seq_len(nrow(settings))) {
        k <- k + 1L
        out[[k]] <- cbind(
            data.frame(image_path = records$image_path[i], setting = settings$name[j],
                       stringsAsFactors = FALSE),
            measure(profiles[[i]], settings$name[j], settings$spar[j]))
    }
    if (i %% 250 == 0) { cat(sprintf("%d / %d\n", i, length(profiles))); flush.console() }
}
out <- do.call(rbind, out)
write.csv(out, file.path(output_dir, "smoothing_variants_per_image.csv"), row.names = FALSE, na = "")

base <- out[out$setting == "spar_0.10", ]
stopifnot(max(abs(base$period_px * 0.0133359375 - records$pitch_1d_um)) < 1e-9)
cat("baseline setting reproduces profiles_per_image.csv\n")
print(aggregate(df ~ setting, data = out, FUN = function(v) round(median(v), 1)))
