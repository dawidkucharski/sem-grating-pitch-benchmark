#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(imager))
suppressPackageStartupMessages(library(pracma))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE)
revision_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
project_dir <- dirname(revision_dir)
output_dir <- file.path(revision_dir, "analysis")
metadata <- read.csv(file.path(output_dir, "sem_metadata.csv"), stringsAsFactors = FALSE)
smoke <- "--smoke" %in% commandArgs(trailingOnly = TRUE)
if (smoke) metadata <- metadata[metadata$detector == "SED", ][1:2, ]

extract_profile <- function(path) {
    image <- load.image(path)
    if (dim(image)[4] >= 3) {
        stopifnot(max(abs(R(image) - G(image))) < 1e-12,
                  max(abs(R(image) - B(image))) < 1e-12)
    }
    matrix_image <- as.matrix(image[, , 1])
    rowMeans(matrix_image)
}

measure_profile <- function(profile, scale, spar = 0.1) {
    fit <- if (is.na(spar)) {
        smooth.spline(seq_along(profile), profile, cv = FALSE)
    } else {
        smooth.spline(seq_along(profile), profile, spar = spar)
    }
    peaks <- findpeaks(-fit$y, nups = 1, ndowns = 1, threshold = 0.02)
    positions <- if (is.null(peaks)) numeric() else sort(peaks[, 2])
    spacings <- diff(positions) * scale
    pitch <- if (length(spacings)) mean(spacings) else NA_real_
    list(
        pitch = pitch,
        positions = positions,
        spacings = spacings,
        min_ratio = if (length(spacings)) min(spacings) / pitch else NA_real_,
        max_ratio = if (length(spacings)) max(spacings) / pitch else NA_real_,
        smoothed = fit$y,
        spar = fit$spar
    )
}

accept <- function(result, lower = 0.6, upper = Inf) {
    is.finite(result$pitch) && length(result$spacings) >= 1 &&
        result$min_ratio >= lower && result$max_ratio <= upper
}

settings <- data.frame(
    name = c("spar_0.01", "spar_0.05", "baseline", "spar_0.2", "spar_0.4", "spar_0.7", "GCV"),
    spar = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.7, NA_real_)
)
records <- vector("list", nrow(metadata))
variants <- vector("list", nrow(metadata) * nrow(settings))
profiles <- vector("list", nrow(metadata))
smoothed_profiles <- vector("list", nrow(metadata))
variant_index <- 0L
for (image_index in seq_len(nrow(metadata))) {
    entry <- metadata[image_index, ]
    profile <- extract_profile(file.path(project_dir, entry$image_path))
    stopifnot(length(profile) == entry$width_px)
    profiles[[image_index]] <- profile
    for (setting_index in seq_len(nrow(settings))) {
        result <- measure_profile(profile, entry$scale_x_um_per_px, settings$spar[setting_index])
        variant_index <- variant_index + 1L
        variants[[variant_index]] <- data.frame(
            image_index = image_index,
            group_label = entry$group_label,
            detector = entry$detector,
            setting = settings$name[setting_index],
            selected_spar = result$spar,
            pitch_um = result$pitch,
            accepted = accept(result)
        )
        if (settings$name[setting_index] == "baseline") {
            tile <- regmatches(entry$image_path, regexec("X([0-9]+)_Y([0-9]+)", entry$image_path))[[1]]
            records[[image_index]] <- data.frame(
                image_index = image_index,
                group_label = entry$group_label,
                detector = entry$detector,
                image_path = entry$image_path,
                tile_x = as.integer(tile[2]),
                tile_y = as.integer(tile[3]),
                pitch_1d_um = result$pitch,
                n_spacings = length(result$spacings),
                sd_spacings_um = if (length(result$spacings) > 1) sd(result$spacings) else NA_real_,
                min_spacing_ratio = result$min_ratio,
                max_spacing_ratio = result$max_ratio,
                accepted = accept(result),
                reason = if (!is.finite(result$pitch)) "fewer_than_two_minima" else if (!accept(result)) "spacing_below_060_mean" else "accepted",
                minima_pixels_1based = paste(result$positions, collapse = ";"),
                stage_x_um = entry$stage_x_um,
                stage_y_um = entry$stage_y_um,
                timestamp = entry$timestamp
            )
            smoothed_profiles[[image_index]] <- result$smoothed
        }
    }
    if (image_index %% 100 == 0) {
        cat(sprintf("Profiles and seven smoothing settings: %d / %d\n", image_index, nrow(metadata)))
        flush.console()
    }
}
records <- do.call(rbind, records)
variants <- do.call(rbind, variants)

if (smoke) {
    old <- read.csv(file.path(project_dir, "R_output_dir", "per_image_pitch_mag400_SED.csv"))
    stopifnot(all(records$accepted))
    for (image_index in seq_len(nrow(records))) {
        expected <- old$pitch_1d_um[basename(old$image_path) == basename(records$image_path[image_index])]
        stopifnot(length(expected) == 1, abs(expected - records$pitch_1d_um[image_index]) < 1e-6)
    }
    print(records[, c("image_path", "pitch_1d_um", "n_spacings", "accepted")], row.names = FALSE)
    cat("SMOKE PASS: horizontal profile has 1280 samples and reproduces the archived per-image pitches.\n")
    quit(status = 0)
}

write.csv(records, file.path(output_dir, "profiles_per_image.csv"), row.names = FALSE, na = "")
write.csv(variants, file.path(output_dir, "smoothing_per_image.csv"), row.names = FALSE, na = "")
saveRDS(list(raw = profiles, smoothed = smoothed_profiles), file.path(output_dir, "profiles.rds"))

sensitivity <- list()
for (group in sort(unique(records$group_label))) {
    for (detector in c("SED", "BED-S")) {
        baseline <- records[records$group_label == group & records$detector == detector, ]
        baseline_pitch <- median(baseline$pitch_1d_um[baseline$accepted])
        for (setting in settings$name) {
            subset <- variants[variants$group_label == group & variants$detector == detector & variants$setting == setting, ]
            common <- merge(subset[subset$accepted, ], baseline[baseline$accepted, c("image_index", "pitch_1d_um")], by = "image_index")
            sensitivity[[length(sensitivity) + 1L]] <- data.frame(
                group_label = group, detector = detector, setting = setting,
                acquired = nrow(subset), accepted = sum(subset$accepted),
                median_um = median(subset$pitch_um[subset$accepted]),
                median_change_um = median(subset$pitch_um[subset$accepted]) - baseline_pitch,
                common_images = nrow(common),
                common_median_absolute_change_um = median(abs(common$pitch_um - common$pitch_1d_um)),
                common_p95_absolute_change_um = unname(quantile(abs(common$pitch_um - common$pitch_1d_um), 0.95))
            )
        }
    }
}
write.csv(do.call(rbind, sensitivity), file.path(output_dir, "smoothing_summary.csv"), row.names = FALSE, na = "")

rejection <- list()
for (group in sort(unique(records$group_label))) {
    for (detector in c("SED", "BED-S")) {
        subset <- records[records$group_label == group & records$detector == detector, ]
        for (lower in c(0, 0.4, 0.5, 0.6, 0.7, 0.8)) {
            for (upper in if (lower == 0.6) c(Inf, 1.4, 1.6, 1.8) else Inf) {
                keep <- is.finite(subset$pitch_1d_um) & subset$min_spacing_ratio >= lower & subset$max_spacing_ratio <= upper
                rejection[[length(rejection) + 1L]] <- data.frame(
                    group_label = group, detector = detector, lower = lower,
                    upper = if (is.infinite(upper)) "none" else as.character(upper),
                    acquired = nrow(subset), accepted = sum(keep),
                    median_um = median(subset$pitch_1d_um[keep]),
                    mad_um = mad(subset$pitch_1d_um[keep])
                )
            }
        }
    }
}
write.csv(do.call(rbind, rejection), file.path(output_dir, "rejection_sensitivity.csv"), row.names = FALSE, na = "")

records$deviation_um <- ave(records$pitch_1d_um, records$group_label, records$detector, FUN = function(values) abs(values - median(values, na.rm = TRUE)))
examples <- unlist(lapply(split(records[records$accepted, ], records$group_label[records$accepted]), function(subset) {
    subset$image_index[order(subset$deviation_um, decreasing = TRUE)[1:2]]
}))
examples <- unique(c(examples, which(records$accepted & records$deviation_um < 0.002)[1:2], which(!records$accepted)[1:2]))
example_profiles <- do.call(rbind, lapply(examples, function(image_index) data.frame(
    image_index = image_index, pixel_1based = seq_along(profiles[[image_index]]),
    intensity = profiles[[image_index]], smoothed = smoothed_profiles[[image_index]]
)))
write.csv(example_profiles, file.path(output_dir, "example_profiles.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "R_session_info.txt"))
print(with(records, table(group_label, detector, accepted)))
cat("FULL PROFILE AUDIT COMPLETE: no previously generated summary fields were reused.\n")