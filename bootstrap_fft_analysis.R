# =========================================================================
# bootstrap_fft_analysis.R
# Bootstrap confidence intervals + 2D FFT pitch validation
# for grating-based SEM magnification assessment
# =========================================================================
#
# This script adds two advanced analyses to the existing pipeline:
#   1. BOOTSTRAP: Non-parametric BCa confidence intervals on all key
#      quantities (median spacing, expanded uncertainty, magnification
#      discrepancy, Bland-Altman bias and LoA). Replaces the internally
#      inconsistent MAD→normal scaling with fully non-parametric inference.
#   2. 2D FFT: Independent pitch determination via 2D Fourier transform,
#      conceptually linked to optical diffractometry and providing an
#      internal cross-validation of the 1D profile-minima method.
#
# Usage: source this file, then call:
#   results <- run_bootstrap_fft_pipeline()
# This will process all four magnifications (400, 500, 750, 1000) for
# SED and BED-S detectors, generate CSV files, and produce publication-
# quality figures.
# =========================================================================

library(imager)
library(splines)
library(pracma)
library(boot)       # for bootstrap confidence intervals
library(parallel)
library(pbapply)

# -------------------------------------------------------------------
# Configuration (mirrors main analysis)
# -------------------------------------------------------------------
base_folder <- "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"

config <- list(
  scale_per_pixel_um      = 0.01333594,
  min_groove_fraction     = 0.6,
  spline_spar             = 0.1,
  k_coverage              = 3,
  typeB_uncertainty_um    = 0.01333594,
  nominal_pitch_um        = 1000 / 600,       # 1.666... um
  n_bootstrap             = 10000,
  confidence_level        = 0.95
)

output_dir <- file.path(base_folder, "R_output_dir")
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)

# =========================================================================
# PART 1: 2D FFT PITCH DETERMINATION
# =========================================================================
# For each image, the 2D FFT magnitude spectrum is computed. The grating
# periodicity produces a characteristic peak pair at spatial frequencies
# +/- 1/T along the direction perpendicular to the grooves. The peak
# position is converted to physical pitch. The FWHM of the fundamental
# peak quantifies groove-period uniformity (sharper = more uniform).
#
# This is the direct SEM analogue of optical diffraction metrology
# (Chernoff2008, Misumi2017): both methods extract pitch from the
# angular position of a diffraction peak, but the SEM-FFT operates
# on the image-intensity Fourier transform rather than on a physical
# diffraction pattern.

#' Compute 2D FFT magnitude spectrum with DC centered
#'
#' @param gray_mat 2D numeric matrix (grayscale image)
#' @return list with: magnitude (centred spectrum), kx, ky (frequency vectors)
fft2d_spectrum <- function(gray_mat) {
  n_rows <- nrow(gray_mat)
  n_cols <- ncol(gray_mat)

  # Apply 2D Hann window to suppress edge artefacts
  hann_x <- 0.5 * (1 - cos(2 * pi * (1:n_cols) / (n_cols + 1)))
  hann_y <- 0.5 * (1 - cos(2 * pi * (1:n_rows) / (n_rows + 1)))
  window <- outer(hann_y, hann_x)
  mat_win <- gray_mat * window

  # 2D FFT via row-column decomposition
  fft_rows <- t(apply(mat_win, 1, fft))
  fft_2d   <- apply(fft_rows, 2, fft)

  # Magnitude spectrum, DC-centred
  mag <- Mod(fft_2d)
  # fftshift: swap quadrants
  mag_centred <- mag
  half_r <- floor(n_rows / 2)
  half_c <- floor(n_cols / 2)
  # We'll use the raw magnitude and extract peaks from the appropriate quadrant
  # For fftshift-equivalent indexing, we work in the range [1, N]
  # Peak at index k corresponds to spatial frequency (k-1)/N cycles/pixel
  # For k > N/2, the physical frequency is (k-1-N)/N (negative frequency)

  # Frequency axes (cycles/pixel)
  kx <- (0:(n_cols - 1)) / n_cols
  ky <- (0:(n_rows - 1)) / n_rows

  list(magnitude = mag, kx = kx, ky = ky,
       n_rows = n_rows, n_cols = n_cols)
}

#' Extract grating pitch from 2D FFT spectrum
#'
#' The grating grooves run approximately horizontally in the image
#' (the 1D profile is obtained by rowMeans, collapsing x). The pitch
#' is therefore along the y-direction, and the FFT fundamental peak
#' appears at ky > 0, kx = 0 (column 1 of the unshifted FFT).
#'
#' A Hann window is applied before the FFT, but some DC leakage into
#' the very-low-frequency bins persists. The search therefore excludes
#' frequencies below 0.003 cyc/px, which corresponds to features
#' larger than ~4.4 µm — well above the grating pitch of ~1.67 µm.
#'
#' @param spec Output of fft2d_spectrum()
#' @param scale_um_per_px Pixel scale in µm
#' @return list with: pitch_um, peak_freq, peak_magnitude,
#'   peak_width_fwhm, snr, quality_flag
extract_fft_pitch <- function(spec, scale_um_per_px) {
  mag    <- spec$magnitude
  n_rows <- spec$n_rows
  n_cols <- spec$n_cols

  # Mask only the 3x3 DC corner, not entire rows/columns
  mag_masked <- mag
  mag_masked[1:3, 1:3] <- 0

  # Search column 1 (kx = 0) for the grating fundamental.
  # Exclude near-DC bins: rows 1:5 correspond to |freq| < 0.0039 cyc/px
  # (pitch > 3.4 µm), which is > 2x the nominal grating pitch.
  # Also exclude the Nyquist row (no physical grating pitch there).
  ky_col <- mag_masked[, 1]
  min_row <- 6
  max_row <- n_rows - 4  # exclude the very last rows (near-DC leakage at neg freq)

  if (min_row > max_row) {
    return(list(pitch_um = NA_real_, peak_freq = NA_real_,
                peak_magnitude = NA_real_, peak_width_fwhm = NA_real_,
                snr = NA_real_, quality_flag = "no_peak_found"))
  }

  search_range <- min_row:max_row
  ky_search <- ky_col[search_range]

  if (max(ky_search) < 1e-10) {
    return(list(pitch_um = NA_real_, peak_freq = NA_real_,
                peak_magnitude = NA_real_, peak_width_fwhm = NA_real_,
                snr = NA_real_, quality_flag = "no_peak_found"))
  }

  peak_rel_idx <- which.max(ky_search)
  peak_row <- search_range[peak_rel_idx]
  peak_val <- ky_search[peak_rel_idx]

  # Convert row index to spatial frequency (cycles/pixel)
  if (peak_row <= n_rows / 2 + 1) {
    freq_ky <- (peak_row - 1) / n_rows
  } else {
    freq_ky <- (peak_row - 1 - n_rows) / n_rows
  }
  freq_ky <- abs(freq_ky)

  if (freq_ky < 1e-10) {
    return(list(pitch_um = NA_real_, peak_freq = NA_real_,
                peak_magnitude = NA_real_, peak_width_fwhm = NA_real_,
                snr = NA_real_, quality_flag = "no_peak_found"))
  }

  pitch_px <- 1 / freq_ky
  pitch_um <- pitch_px * scale_um_per_px

  # Peak width (FWHM) from the 1D slice at kx=0 (column 1)
  half_max <- peak_val / 2
  above <- which(ky_search >= half_max)
  if (length(above) >= 2) {
    fwhm_py <- (search_range[above[length(above)]] -
                search_range[above[1]]) / n_rows
  } else {
    fwhm_py <- NA_real_
  }

  # Signal-to-noise ratio: peak vs median of search region
  noise_floor <- median(ky_search[ky_search > 0], na.rm = TRUE)
  if (is.na(noise_floor) || noise_floor == 0) noise_floor <- 1
  snr <- peak_val / noise_floor

  quality <- if (snr > 10) "excellent" else if (snr > 5) "good" else "marginal"

  list(pitch_um         = pitch_um,
       peak_freq        = freq_ky,
       peak_magnitude   = as.numeric(peak_val),
       peak_width_fwhm  = fwhm_py,
       snr              = as.numeric(snr),
       quality_flag     = quality)
}

# Helper for NULL-safe fallback (must be defined before use)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Process a single image: 1D profile pitch + 2D FFT pitch
#'
#' @param image_path Full path to image file
#' @param global_mean_intensity Mean intensity for normalisation
#' @return data.frame with per-image metrics (1 row)
process_image_dual <- function(image_path, global_mean_intensity) {
  # Load image
  img <- tryCatch(load.image(image_path), error = function(e) NULL)
  if (is.null(img)) return(NULL)

  # Grayscale conversion
  if (dim(img)[3] > 1) {
    gray_img <- grayscale(img)
  } else {
    gray_img <- img
  }

  # Intensity normalisation
  gray_img <- gray_img + (global_mean_intensity - mean(gray_img))
  gray_mat <- as.matrix(gray_img[, , 1])

  # ---- 1D profile method (existing) ----
  intensity_profile <- rowMeans(gray_mat)
  spline_fit <- smooth.spline(
    x = 1:length(intensity_profile),
    y = intensity_profile,
    spar = config$spline_spar
  )
  smoothed_profile <- spline_fit$y

  # Detect minima
  inverted_profile <- -smoothed_profile
  minima <- tryCatch(
    findpeaks(inverted_profile, nups = 1, ndowns = 1, threshold = 0.02),
    error = function(e) NULL
  )
  if (is.null(minima)) return(NULL)

  minima_pos <- minima[, 2]
  if (length(minima_pos) < 2) return(NULL)

  dist_px <- diff(minima_pos)
  dist_um <- dist_px * config$scale_per_pixel_um

  # Rejection criterion
  mean_dist <- mean(dist_um, na.rm = TRUE)
  if (any(dist_um < mean_dist * config$min_groove_fraction)) return(NULL)

  pitch_1d_um <- mean(dist_um, na.rm = TRUE)
  n_grooves   <- length(dist_um)
  sd_grooves_um <- sd(dist_um, na.rm = TRUE)

  # Profile height (min-to-max of smoothed profile)
  profile_height <- max(smoothed_profile) - min(smoothed_profile)

  # Detector tag
  fname <- basename(image_path)
  detector <- if (grepl("_SED_", fname, ignore.case = TRUE)) "SED" else "BED-S"

  # Tile indices from filename
  m <- regexpr("X[0-9]+_Y[0-9]+", fname)
  tile_x <- NA_integer_
  tile_y <- NA_integer_
  if (m != -1) {
    xy_str <- regmatches(fname, m)
    parts  <- strsplit(xy_str, "_")[[1]]
    tile_x <- as.integer(sub("X", "", parts[1]))
    tile_y <- as.integer(sub("Y", "", parts[2]))
  }

  # ---- 2D FFT method (new) ----
  fft_spec <- fft2d_spectrum(gray_mat)
  fft_result <- extract_fft_pitch(fft_spec, config$scale_per_pixel_um)

  # ---- Magnification tag ----
  mag_tag <- NA_integer_
  for (m in c(400, 500, 750, 1000)) {
    if (grepl(sprintf("mag%d", m), image_path)) mag_tag <- m
  }

  data.frame(
    magnification       = mag_tag,
    detector            = detector,
    image_path          = image_path,
    tile_x              = tile_x,
    tile_y              = tile_y,
    pitch_1d_um         = round(pitch_1d_um, 6),
    n_grooves           = n_grooves,
    sd_grooves_um       = round(sd_grooves_um, 6),
    profile_height      = round(profile_height, 6),
    pitch_fft_um        = round(fft_result$pitch_um, 6),
    fft_peak_freq       = round(fft_result$peak_freq, 8),
    fft_snr             = round(fft_result$snr, 2),
    fft_quality         = fft_result$quality_flag,
    fft_peak_width_fwhm = round(fft_result$peak_width_fwhm %||% NA_real_, 8),
    stringsAsFactors    = FALSE
  )
}

# =========================================================================
# PART 2: BOOTSTRAP INFERENCE
# =========================================================================

#' Bootstrap confidence intervals for key metrological quantities
#'
#' Uses BCa (bias-corrected and accelerated) bootstrap on per-image
#' spacing values. All inference is non-parametric — no distributional
#' assumptions.
#'
#' @param x Numeric vector of per-image mean spacing values (µm)
#' @param typeB Type B uncertainty (µm)
#' @param k Coverage factor
#' @param d_nom Nominal pitch (µm)
#' @param R Number of bootstrap replicates
#' @param conf Confidence level
#' @return list with point estimates and BCa confidence intervals
bootstrap_metrology <- function(x, typeB, k, d_nom, R = 10000, conf = 0.95) {
  n <- length(x)
  if (n < 5) return(NULL)

  # Statistic function for boot(): returns vector of derived quantities
  stat_fun <- function(data, indices) {
    d <- data[indices]
    med     <- median(d)
    mad_val <- mad(d, constant = 1.4826)
    uA      <- mad_val / sqrt(length(d))
    uC      <- sqrt(uA^2 + typeB^2)
    U       <- k * uC
    eps_M   <- (med - d_nom) / d_nom * 100
    c(median_um = med, uA_um = uA, uC_um = uC, U_um = U, eps_M_pct = eps_M)
  }

  boot_res <- boot(data = x, statistic = stat_fun, R = R, parallel = "multicore",
                   ncpus = max(1, detectCores() - 1))

  # Extract BCa CIs for each quantity
  quantities <- c("median_um", "uA_um", "uC_um", "U_um", "eps_M_pct")
  ci_list <- lapply(quantities, function(q) {
    idx <- which(q == quantities)
    tryCatch({
      ci <- boot.ci(boot_res, conf = conf, type = "bca", index = idx)
      if (!is.null(ci) && !is.null(ci$bca)) {
        c(lower = ci$bca[4], upper = ci$bca[5])
      } else {
        # fallback to percentile
        c(lower = quantile(boot_res$t[, idx], (1 - conf) / 2, na.rm = TRUE),
          upper = quantile(boot_res$t[, idx], 1 - (1 - conf) / 2, na.rm = TRUE))
      }
    }, error = function(e) {
      c(lower = quantile(boot_res$t[, idx], (1 - conf) / 2, na.rm = TRUE),
        upper = quantile(boot_res$t[, idx], 1 - (1 - conf) / 2, na.rm = TRUE))
    })
  })
  names(ci_list) <- quantities

  # Point estimates
  point <- stat_fun(x, 1:n)

  list(
    point_estimates = list(
      median_um = point[1],
      uA_um     = point[2],
      uC_um     = point[3],
      U_um      = point[4],
      eps_M_pct = point[5]
    ),
    ci = ci_list,
    n_images = n,
    boot_obj = boot_res
  )
}

#' Bootstrap Bland-Altman statistics
#'
#' @param sed Vector of SED spacing values
#' @param beds Vector of BED-S spacing values (paired, same length as sed)
#' @param R Number of bootstrap replicates
#' @param conf Confidence level
bootstrap_bland_altman <- function(sed, beds, R = 10000, conf = 0.95) {
  n <- length(sed)
  if (n < 5) return(NULL)

  stat_fun <- function(data, indices) {
    s <- data[indices, 1]
    b <- data[indices, 2]
    diff <- s - b
    bias  <- mean(diff)
    sd_diff <- sd(diff)
    loa_lower <- bias - 1.96 * sd_diff
    loa_upper <- bias + 1.96 * sd_diff
    c(bias_um = bias, loa_lower_um = loa_lower, loa_upper_um = loa_upper,
      sd_diff_um = sd_diff)
  }

  paired_data <- cbind(sed, beds)
  boot_res <- boot(data = paired_data, statistic = stat_fun, R = R,
                   parallel = "multicore", ncpus = max(1, detectCores() - 1))

  quantities <- c("bias_um", "loa_lower_um", "loa_upper_um", "sd_diff_um")
  ci_list <- lapply(seq_along(quantities), function(idx) {
    tryCatch({
      ci <- boot.ci(boot_res, conf = conf, type = "bca", index = idx)
      if (!is.null(ci) && !is.null(ci$bca)) {
        c(lower = ci$bca[4], upper = ci$bca[5])
      } else {
        c(lower = quantile(boot_res$t[, idx], (1 - conf) / 2, na.rm = TRUE),
          upper = quantile(boot_res$t[, idx], 1 - (1 - conf) / 2, na.rm = TRUE))
      }
    }, error = function(e) {
      c(lower = quantile(boot_res$t[, idx], (1 - conf) / 2, na.rm = TRUE),
        upper = quantile(boot_res$t[, idx], 1 - (1 - conf) / 2, na.rm = TRUE))
    })
  })
  names(ci_list) <- quantities

  point <- stat_fun(paired_data, 1:n)
  list(
    point_estimates = list(
      bias_um      = point[1],
      loa_lower_um = point[2],
      loa_upper_um = point[3],
      sd_diff_um   = point[4]
    ),
    ci = ci_list,
    n_pairs = n,
    boot_obj = boot_res
  )
}

# =========================================================================
# PART 3: MAIN PIPELINE
# =========================================================================

#' Run full bootstrap + FFT analysis pipeline
#'
#' Processes images for all four magnifications and both detectors,
#' computes per-image 1D and FFT pitch values, then runs bootstrap
#' inference on all key quantities.
#'
#' @param mags Vector of magnifications to process
#' @param max_images_per_combo Max images per mag-detector combo (NULL = all)
#' @return Invisible list of results; also writes CSV and PDF files
run_bootstrap_fft_pipeline <- function(mags = c(400, 500, 750, 1000),
                                        max_images_per_combo = NULL) {
  all_per_image <- list()
  bootstrap_summaries <- list()
  bland_altman_summaries <- list()

  for (mag in mags) {
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat(sprintf("Processing magnification %d\n", mag))
    cat(paste(rep("=", 70), collapse = ""), "\n")

    for (det in c("SED", "BED-S")) {
      cat(sprintf("\n--- %s detector ---\n", det))

      # Build folder path
      if (det == "SED") {
        img_folder <- file.path(base_folder,
                                sprintf("mag%d_SED-BED", mag), "SED")
      } else {
        img_folder <- file.path(base_folder,
                                sprintf("mag%d_SED-BED", mag), "BED-S")
      }

      img_files <- list.files(img_folder, full.names = TRUE,
                              pattern = "\\.(jpg|png|tif|bmp)$")
      if (length(img_files) == 0) {
        cat(sprintf("  No images found in %s\n", img_folder))
        next
      }

      # Optionally subsample
      if (!is.null(max_images_per_combo) &&
          length(img_files) > max_images_per_combo) {
        set.seed(123)
        img_files <- sample(img_files, max_images_per_combo)
      }

      cat(sprintf("  Processing %d images...\n", length(img_files)))

      # Compute global mean intensity for normalisation
      cl <- makeCluster(max(1, detectCores() - 1))
      clusterExport(cl, c("config", "fft2d_spectrum", "extract_fft_pitch",
                          "process_image_dual"),
                    envir = environment())
      clusterEvalQ(cl, {
        library(imager)
        library(splines)
        library(pracma)
      })

      mean_intensities <- parSapply(cl, img_files, function(fp) {
        img <- tryCatch(load.image(fp), error = function(e) NULL)
        if (is.null(img)) return(NA_real_)
        if (dim(img)[3] > 1) img <- grayscale(img)
        mean(as.matrix(img[, , 1]))
      })
      stopCluster(cl)

      global_mean <- mean(mean_intensities, na.rm = TRUE)
      cat(sprintf("  Global mean intensity: %.3f\n", global_mean))

      # Process images (sequential with pbapply for progress)
      results_list <- pblapply(img_files, process_image_dual,
                               global_mean_intensity = global_mean)
      results_list <- results_list[!sapply(results_list, is.null)]

      if (length(results_list) == 0) {
        cat("  No valid images processed.\n")
        next
      }

      per_image_df <- do.call(rbind, results_list)
      per_image_df$magnification <- mag
      per_image_df$detector <- det

      combo_key <- sprintf("mag%d_%s", mag, det)
      all_per_image[[combo_key]] <- per_image_df

      # Save per-image data
      per_img_file <- file.path(output_dir,
                                sprintf("per_image_pitch_mag%d_%s.csv", mag, det))
      write.csv(per_image_df, per_img_file, row.names = FALSE)
      cat(sprintf("  Per-image data saved: %s (%d rows)\n",
                  per_img_file, nrow(per_image_df)))

      # ---- Bootstrap inference ----
      cat(sprintf("  Running bootstrap (R=%d)...\n", config$n_bootstrap))
      pitch_values <- per_image_df$pitch_1d_um

      boot_results <- bootstrap_metrology(
        x     = pitch_values,
        typeB = config$typeB_uncertainty_um,
        k     = config$k_coverage,
        d_nom = config$nominal_pitch_um,
        R     = config$n_bootstrap,
        conf  = config$confidence_level
      )

      if (!is.null(boot_results)) {
        bootstrap_summaries[[combo_key]] <- boot_results

        # Print results
        cat(sprintf("  Median spacing: %.3f µm [%.3f, %.3f] (95%% BCa)\n",
                    boot_results$point_estimates$median_um,
                    boot_results$ci$median_um[1],
                    boot_results$ci$median_um[2]))
        cat(sprintf("  Expanded U (k=%d): %.3f µm [%.3f, %.3f]\n",
                    config$k_coverage,
                    boot_results$point_estimates$U_um,
                    boot_results$ci$U_um[1],
                    boot_results$ci$U_um[2]))
        cat(sprintf("  ε_M: %.2f%% [%.2f, %.2f]\n",
                    boot_results$point_estimates$eps_M_pct,
                    boot_results$ci$eps_M_pct[1],
                    boot_results$ci$eps_M_pct[2]))
      }

      cat(sprintf("  FFT pitch (median): %.3f µm (mean: %.3f µm)\n",
                  median(per_image_df$pitch_fft_um, na.rm = TRUE),
                  mean(per_image_df$pitch_fft_um, na.rm = TRUE)))
      cat(sprintf("  1D-FFT agreement (median diff): %.4f µm\n",
                  median(per_image_df$pitch_1d_um - per_image_df$pitch_fft_um,
                         na.rm = TRUE)))
      cat(sprintf("  FFT quality: %s\n",
                  paste(sprintf("%s=%d", names(table(per_image_df$fft_quality)),
                                table(per_image_df$fft_quality)),
                        collapse = ", ")))
    }

    # ---- Bland-Altman bootstrap at this magnification ----
    sed_key  <- sprintf("mag%d_SED", mag)
    beds_key <- sprintf("mag%d_BED-S", mag)

    if (!is.null(all_per_image[[sed_key]]) &&
        !is.null(all_per_image[[beds_key]])) {
      sed_df  <- all_per_image[[sed_key]]
      beds_df <- all_per_image[[beds_key]]

      # Pair by tile position
      paired <- merge(sed_df[, c("tile_x", "tile_y", "pitch_1d_um")],
                      beds_df[, c("tile_x", "tile_y", "pitch_1d_um")],
                      by = c("tile_x", "tile_y"),
                      suffixes = c("_sed", "_beds"))
      paired <- paired[complete.cases(paired), ]

      if (nrow(paired) >= 5) {
        cat(sprintf("\n  Bland-Altman bootstrap: %d paired tiles\n", nrow(paired)))
        ba_boot <- bootstrap_bland_altman(
          sed  = paired$pitch_1d_um_sed,
          beds = paired$pitch_1d_um_beds,
          R    = config$n_bootstrap,
          conf = config$confidence_level
        )

        if (!is.null(ba_boot)) {
          ba_key <- sprintf("mag%d", mag)
          bland_altman_summaries[[ba_key]] <- ba_boot

          cat(sprintf("  Bias: %.4f µm [%.4f, %.4f]\n",
                      ba_boot$point_estimates$bias_um,
                      ba_boot$ci$bias_um[1],
                      ba_boot$ci$bias_um[2]))
          cat(sprintf("  LoA: [%.4f, %.4f] µm\n",
                      ba_boot$point_estimates$loa_lower_um,
                      ba_boot$point_estimates$loa_upper_um))
        }
      } else {
        cat(sprintf("\n  Insufficient paired tiles for BA bootstrap (%d)\n",
                    nrow(paired)))
      }
    }
  }

  # =========================================================================
  # Export bootstrap summary tables
  # =========================================================================
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("Exporting summary tables...\n")

  # Bootstrap summary for all combos
  boot_summary_rows <- list()
  for (key in names(bootstrap_summaries)) {
    parts <- strsplit(key, "_")[[1]]
    br <- bootstrap_summaries[[key]]
    boot_summary_rows[[length(boot_summary_rows) + 1]] <- data.frame(
      magnification = as.integer(sub("mag", "", parts[1])),
      detector      = parts[2],
      n_images      = br$n_images,
      median_um     = br$point_estimates$median_um,
      median_ci_low = br$ci$median_um[1],
      median_ci_upp = br$ci$median_um[2],
      uA_um         = br$point_estimates$uA_um,
      uA_ci_low     = br$ci$uA_um[1],
      uA_ci_upp     = br$ci$uA_um[2],
      uC_um         = br$point_estimates$uC_um,
      uC_ci_low     = br$ci$uC_um[1],
      uC_ci_upp     = br$ci$uC_um[2],
      U_um          = br$point_estimates$U_um,
      U_ci_low      = br$ci$U_um[1],
      U_ci_upp      = br$ci$U_um[2],
      eps_M_pct     = br$point_estimates$eps_M_pct,
      eps_M_ci_low  = br$ci$eps_M_pct[1],
      eps_M_ci_upp  = br$ci$eps_M_pct[2],
      stringsAsFactors = FALSE
    )
  }
  boot_summary_df <- do.call(rbind, boot_summary_rows)
  boot_file <- file.path(output_dir, "bootstrap_summary.csv")
  write.csv(boot_summary_df, boot_file, row.names = FALSE)
  cat(sprintf("Bootstrap summary: %s\n", boot_file))

  # Bland-Altman bootstrap summary
  ba_summary_rows <- list()
  for (key in names(bland_altman_summaries)) {
    ba <- bland_altman_summaries[[key]]
    ba_summary_rows[[length(ba_summary_rows) + 1]] <- data.frame(
      magnification   = as.integer(sub("mag", "", key)),
      n_pairs         = ba$n_pairs,
      bias_um         = ba$point_estimates$bias_um,
      bias_ci_low     = ba$ci$bias_um[1],
      bias_ci_upp     = ba$ci$bias_um[2],
      loa_lower_um    = ba$point_estimates$loa_lower_um,
      loa_lower_ci_low = ba$ci$loa_lower_um[1],
      loa_lower_ci_upp = ba$ci$loa_lower_um[2],
      loa_upper_um    = ba$point_estimates$loa_upper_um,
      loa_upper_ci_low = ba$ci$loa_upper_um[1],
      loa_upper_ci_upp = ba$ci$loa_upper_um[2],
      sd_diff_um      = ba$point_estimates$sd_diff_um,
      sd_diff_ci_low  = ba$ci$sd_diff_um[1],
      sd_diff_ci_upp  = ba$ci$sd_diff_um[2],
      stringsAsFactors = FALSE
    )
  }
  ba_summary_df <- do.call(rbind, ba_summary_rows)
  ba_file <- file.path(output_dir, "bland_altman_bootstrap_summary.csv")
  write.csv(ba_summary_df, ba_file, row.names = FALSE)
  cat(sprintf("Bland-Altman bootstrap summary: %s\n", ba_file))

  # =========================================================================
  # Generate publication-quality figures
  # =========================================================================
  cat("\nGenerating figures...\n")

  # Figure 1: 1D vs FFT pitch agreement (Bland-Altman style)
  cat("  Figure: 1D-FFT pitch agreement...\n")
  all_combined <- do.call(rbind, all_per_image)
  # Remove rows where FFT failed
  valid_fft <- all_combined[!is.na(all_combined$pitch_fft_um) &
                            !is.na(all_combined$pitch_1d_um), ]
  if (nrow(valid_fft) > 10) {
    diff_fft <- valid_fft$pitch_1d_um - valid_fft$pitch_fft_um
    mean_diff <- mean(diff_fft)
    sd_diff   <- sd(diff_fft)

    pdf(file.path(base_folder, "fft_1d_agreement.pdf"), width = 8, height = 6)
    par(mar = c(5, 5, 4, 2))
    avg_pitch <- (valid_fft$pitch_1d_um + valid_fft$pitch_fft_um) / 2
    plot(avg_pitch, diff_fft,
         xlab = "Mean Pitch (1D + FFT) / 2 [µm]",
         ylab = "Difference (1D - FFT) [µm]",
         main = "Agreement: 1D Profile vs 2D FFT Pitch",
         pch = 16, col = rgb(0, 0, 1, 0.3),
         cex.lab = 1.2, cex.main = 1.3)
    abline(h = mean_diff, col = "red", lwd = 2)
    abline(h = mean_diff + 1.96 * sd_diff, col = "blue", lwd = 1.5, lty = 2)
    abline(h = mean_diff - 1.96 * sd_diff, col = "blue", lwd = 1.5, lty = 2)
    abline(h = 0, col = "grey50", lty = 3)
    legend("topright",
           legend = c(sprintf("Mean diff: %.4f µm", mean_diff),
                      sprintf("95%% LoA: [%.4f, %.4f] µm",
                              mean_diff - 1.96 * sd_diff,
                              mean_diff + 1.96 * sd_diff)),
           bty = "n", cex = 0.9)
    dev.off()
    cat(sprintf("    Mean 1D-FFT difference: %.4f µm (SD: %.4f µm)\n",
                mean_diff, sd_diff))
  }

  # Figure 2: Bootstrap CI forest plot for ε_M
  cat("  Figure: Bootstrap ε_M forest plot...\n")
  if (nrow(boot_summary_df) > 0) {
    pdf(file.path(base_folder, "bootstrap_epsM_forest.pdf"), width = 10, height = 5)
    par(mar = c(5, 10, 4, 2))
    n_combos <- nrow(boot_summary_df)
    y_pos <- rev(1:n_combos)
    labels <- sprintf("%dx %s", boot_summary_df$magnification,
                      boot_summary_df$detector)
    
    plot(boot_summary_df$eps_M_pct, y_pos,
         xlim = range(c(boot_summary_df$eps_M_ci_low,
                        boot_summary_df$eps_M_ci_upp), na.rm = TRUE) +
                c(-0.2, 0.2),
         yaxt = "n", ylab = "", xlab = expression(epsilon[M] * " [%]"),
         main = expression("Apparent Magnification Discrepancy " * epsilon[M] *
                          " -- Bootstrap 95% BCa CIs"),
         pch = 16, cex = 1.5, cex.lab = 1.2, cex.main = 1.1)
    arrows(boot_summary_df$eps_M_ci_low, y_pos,
           boot_summary_df$eps_M_ci_upp, y_pos,
           angle = 90, code = 3, length = 0.05, lwd = 2)
    axis(2, at = y_pos, labels = labels, las = 2, cex.axis = 0.9)
    abline(v = 0, lty = 3, col = "grey50")
    dev.off()
  }

  # Figure 3: FFT spectrum example (first valid image from mag=1000)
  cat("  Figure: FFT spectrum example...\n")
  example_img_path <- NULL
  for (mag in c(1000, 750, 500, 400)) {
    combo <- sprintf("mag%d_SED", mag)
    if (!is.null(all_per_image[[combo]]) &&
        nrow(all_per_image[[combo]]) > 0) {
      example_img_path <- all_per_image[[combo]]$image_path[1]
      break
    }
  }

  if (!is.null(example_img_path) && file.exists(example_img_path)) {
    img <- load.image(example_img_path)
    if (dim(img)[3] > 1) img <- grayscale(img)
    gray_mat <- as.matrix(img[, , 1])

    spec <- fft2d_spectrum(gray_mat)

    # Log-magnitude for visualisation
    log_mag <- log1p(spec$magnitude)

    pdf(file.path(base_folder, "fft_spectrum_example.pdf"), width = 12, height = 5)
    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

    # Original image (x = rows = y-direction, y = cols = x-direction)
    image(1:nrow(gray_mat), 1:ncol(gray_mat), gray_mat,
          col = grey(seq(0, 1, length.out = 256)),
          xlab = "Y [px]", ylab = "X [px]",
          main = "Original SEM Image", useRaster = TRUE)

    # FFT magnitude (log scale). DC is at (1,1).
    image(1:spec$n_rows, 1:spec$n_cols,
          log_mag,
          col = hcl.colors(256, "YlOrRd", rev = TRUE),
          xlab = "ky index", ylab = "kx index",
          main = "2D FFT (log magnitude)", useRaster = TRUE)

    # 1D slice along ky at kx = 1 (DC row in unshifted FFT = column 1)
    # This is where the grating fundamental peak appears.
    ky_slice <- spec$magnitude[, 1]
    ky_slice[1:3] <- NA  # mask DC
    plot(spec$ky, ky_slice, type = "l", lwd = 1.5,
         xlab = "Spatial Frequency ky [cycles/px]",
         ylab = "|FFT| Magnitude",
         main = "FFT Slice at kx = 0",
         cex.lab = 1.1)
    # Mark detected peak
    fft_res <- extract_fft_pitch(spec, config$scale_per_pixel_um)
    abline(v = fft_res$peak_freq, col = "red", lwd = 2, lty = 2)
    legend("topright",
           legend = c(sprintf("Peak: ky=%.4f cyc/px", fft_res$peak_freq),
                      sprintf("Pitch: %.3f um", fft_res$pitch_um)),
           col = c("red", NA), lty = c(2, NA), lwd = c(2, NA),
           bty = "n", cex = 0.9)

    dev.off()
  }

  cat("\nDone. All output files in:", output_dir, "\n")

  invisible(list(
    per_image     = all_per_image,
    bootstrap     = bootstrap_summaries,
    bland_altman  = bland_altman_summaries,
    config        = config
  ))
}

# =========================================================================
# Run if sourced directly
# =========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Bootstrap + FFT Analysis Pipeline\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("To run: results <- run_bootstrap_fft_pipeline()\n")
cat("For quick test: results <- run_bootstrap_fft_pipeline(max_images_per_combo = 20)\n")
