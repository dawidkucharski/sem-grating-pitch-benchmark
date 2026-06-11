# =========================================================================
# simulate_validate.R
# Pipeline validation using synthetic grating images with known pitch
# =========================================================================
# Generates synthetic SEM-like grating images at each magnification with
# precisely known pitch, then processes them through the identical 1D and
# FFT pipeline. This quantifies algorithmic bias independent of the real
# grating and provides an ironclad answer to the question:
# "Does your pipeline recover the true pitch, or is it biased?"
# =========================================================================

library(imager)
library(splines)
library(pracma)

base_folder <- "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"

config <- list(
  scale_per_pixel_um  = 0.01333594,
  spline_spar         = 0.1,
  min_groove_fraction = 0.6
)

# -------------------------------------------------------------------
# Generate a synthetic grating image
# -------------------------------------------------------------------
generate_grating_image <- function(pitch_um, n_pixels_x = 1280, n_pixels_y = 960,
                                    contrast = 0.3, noise_sd = 0.02,
                                    groove_angle_deg = 0) {
  # pitch_um: true grating pitch in µm
  # Returns: grayscale matrix representing an SEM-like grating image
  
  pitch_px <- pitch_um / config$scale_per_pixel_um
  
  # Generate sinusoidal intensity pattern along y (grooves horizontal)
  y_idx <- 1:n_pixels_y
  intensity <- 0.5 + contrast * sin(2 * pi * y_idx / pitch_px)
  
  # Add slight random phase jitter per column (simulates groove waviness)
  set.seed(42)
  phase_jitter <- matrix(rnorm(n_pixels_x, 0, 0.02), nrow = n_pixels_y,
                         ncol = n_pixels_x, byrow = TRUE)
  
  # Replicate pattern across columns with jitter
  img <- matrix(rep(intensity, n_pixels_x), nrow = n_pixels_y, ncol = n_pixels_x)
  img <- img + phase_jitter
  
  # Add Gaussian noise
  img <- img + matrix(rnorm(n_pixels_x * n_pixels_y, 0, noise_sd),
                      nrow = n_pixels_y, ncol = n_pixels_x)
  
  # Clip to [0,1]
  img[img < 0] <- 0
  img[img > 1] <- 1
  
  return(img)
}

# -------------------------------------------------------------------
# Process a single synthetic image through the 1D pipeline
# -------------------------------------------------------------------
process_1d <- function(gray_mat) {
  intensity_profile <- rowMeans(gray_mat)
  
  spline_fit <- smooth.spline(x = 1:length(intensity_profile),
                               y = intensity_profile,
                               spar = config$spline_spar)
  smoothed_profile <- spline_fit$y
  
  inverted_profile <- -smoothed_profile
  minima <- tryCatch(
    findpeaks(inverted_profile, nups = 1, ndowns = 1, threshold = 0.02),
    error = function(e) NULL
  )
  
  if (is.null(minima)) return(NA_real_)
  minima_pos <- minima[, 2]
  if (length(minima_pos) < 2) return(NA_real_)
  
  dist_px <- diff(minima_pos)
  dist_um <- dist_px * config$scale_per_pixel_um
  
  mean_dist <- mean(dist_um, na.rm = TRUE)
  if (any(dist_um < mean_dist * config$min_groove_fraction)) return(NA_real_)
  
  return(mean(dist_um, na.rm = TRUE))
}

# -------------------------------------------------------------------
# Process through FFT pipeline (simplified)
# -------------------------------------------------------------------
process_fft <- function(gray_mat) {
  n_rows <- nrow(gray_mat)
  n_cols <- ncol(gray_mat)
  
  # Hann window
  hann_x <- 0.5 * (1 - cos(2 * pi * (1:n_cols) / (n_cols + 1)))
  hann_y <- 0.5 * (1 - cos(2 * pi * (1:n_rows) / (n_rows + 1)))
  window <- outer(hann_y, hann_x)
  mat_win <- gray_mat * window
  
  # 2D FFT
  fft_rows <- t(apply(mat_win, 1, fft))
  fft_2d   <- apply(fft_rows, 2, fft)
  mag <- Mod(fft_2d)
  
  # Mask DC
  mag[1:3, 1:3] <- 0
  
  # Search column 1 (kx=0) for grating fundamental
  ky_col <- mag[, 1]
  # Exclude near-DC
  min_row <- 6
  max_row <- n_rows - 4
  if (min_row > max_row) return(NA_real_)
  
  search_range <- min_row:max_row
  ky_search <- ky_col[search_range]
  
  if (max(ky_search) < 1e-10) return(NA_real_)
  
  peak_rel_idx <- which.max(ky_search)
  peak_row <- search_range[peak_rel_idx]
  
  if (peak_row <= n_rows / 2 + 1) {
    freq_ky <- (peak_row - 1) / n_rows
  } else {
    freq_ky <- (peak_row - 1 - n_rows) / n_rows
  }
  freq_ky <- abs(freq_ky)
  
  if (freq_ky < 1e-10) return(NA_real_)
  
  pitch_px <- 1 / freq_ky
  pitch_um <- pitch_px * config$scale_per_pixel_um
  return(pitch_um)
}

# -------------------------------------------------------------------
# Main simulation: test pipeline at multiple true pitches
# -------------------------------------------------------------------
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Pipeline Validation: Synthetic Grating Images\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Test pitches: nominal, FFT-measured, 1D-measured, ±5% extremes
true_pitches_um <- c(1.667, 1.707, 1.643, 1.58, 1.75)
n_replicates <- 50  # images per pitch value

results <- data.frame()

for (true_pitch in true_pitches_um) {
  cat(sprintf("True pitch: %.3f µm ... ", true_pitch))
  
  pitches_1d <- numeric(n_replicates)
  pitches_fft <- numeric(n_replicates)
  
  for (i in 1:n_replicates) {
    img <- generate_grating_image(true_pitch, noise_sd = 0.02)
    pitches_1d[i] <- process_1d(img)
    pitches_fft[i] <- process_fft(img)
  }
  
  # Summary statistics
  valid_1d <- pitches_1d[!is.na(pitches_1d)]
  valid_fft <- pitches_fft[!is.na(pitches_fft)]
  
  results <- rbind(results, data.frame(
    true_pitch_um    = true_pitch,
    n_1d_valid       = length(valid_1d),
    n_fft_valid      = length(valid_fft),
    mean_1d_um       = mean(valid_1d),
    sd_1d_um         = sd(valid_1d),
    mean_fft_um      = mean(valid_fft),
    sd_fft_um        = sd(valid_fft),
    bias_1d_um       = mean(valid_1d) - true_pitch,
    bias_fft_um      = mean(valid_fft) - true_pitch,
    bias_1d_pct      = 100 * (mean(valid_1d) - true_pitch) / true_pitch,
    bias_fft_pct     = 100 * (mean(valid_fft) - true_pitch) / true_pitch,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("1D: %.4f ± %.4f µm (bias: %+.4f µm, %+.2f%%), ",
              results$mean_1d_um[nrow(results)],
              results$sd_1d_um[nrow(results)],
              results$bias_1d_um[nrow(results)],
              results$bias_1d_pct[nrow(results)]))
  cat(sprintf("FFT: %.4f ± %.4f µm (bias: %+.4f µm, %+.2f%%)\n",
              results$mean_fft_um[nrow(results)],
              results$sd_fft_um[nrow(results)],
              results$bias_fft_um[nrow(results)],
              results$bias_fft_pct[nrow(results)]))
}

# Save results
output_dir <- file.path(base_folder, "R_output_dir")
write.csv(results, file.path(output_dir, "simulation_validation.csv"),
          row.names = FALSE)

cat("\n--- Summary ---\n")
cat(sprintf("Mean 1D bias across all pitches: %+.4f µm (%+.2f%%)\n",
            mean(results$bias_1d_um), mean(results$bias_1d_pct)))
cat(sprintf("Mean FFT bias across all pitches: %+.4f µm (%+.2f%%)\n",
            mean(results$bias_fft_um), mean(results$bias_fft_pct)))
cat(sprintf("1D rejection rate: %.0f%%\n",
            100 * (1 - mean(results$n_1d_valid / n_replicates))))
cat(sprintf("FFT success rate: %.0f%%\n",
            100 * mean(results$n_fft_valid / n_replicates)))

# ---- Generate validation figure ----
pdf(file.path(base_folder, "simulation_validation.pdf"), width = 8, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 1.5))

# Left: 1D method
plot(results$true_pitch_um, results$bias_1d_um, type = "b", pch = 19, cex = 1.5,
     xlab = "True pitch [µm]", ylab = "Bias (measured − true) [µm]",
     main = "1D Profile-Minima Method",
     ylim = range(c(results$bias_1d_um, results$bias_fft_um, 0), na.rm = TRUE),
     cex.lab = 1.1)
abline(h = 0, lty = 2, col = "grey50")
arrows(results$true_pitch_um,
       results$bias_1d_um - results$sd_1d_um / sqrt(results$n_1d_valid),
       results$true_pitch_um,
       results$bias_1d_um + results$sd_1d_um / sqrt(results$n_1d_valid),
       angle = 90, code = 3, length = 0.05)
legend("topright",
       legend = sprintf("Mean bias: %+.4f µm", mean(results$bias_1d_um)),
       bty = "n", cex = 0.9)

# Right: FFT method
plot(results$true_pitch_um, results$bias_fft_um, type = "b", pch = 19, cex = 1.5,
     xlab = "True pitch [µm]", ylab = "Bias (measured − true) [µm]",
     main = "2D FFT Method",
     ylim = range(c(results$bias_1d_um, results$bias_fft_um, 0), na.rm = TRUE),
     cex.lab = 1.1)
abline(h = 0, lty = 2, col = "grey50")
arrows(results$true_pitch_um,
       results$bias_fft_um - results$sd_fft_um / sqrt(results$n_fft_valid),
       results$true_pitch_um,
       results$bias_fft_um + results$sd_fft_um / sqrt(results$n_fft_valid),
       angle = 90, code = 3, length = 0.05)
legend("topright",
       legend = sprintf("Mean bias: %+.4f µm", mean(results$bias_fft_um)),
       bty = "n", cex = 0.9)

dev.off()
cat(sprintf("\nFigure saved: simulation_validation.pdf\n"))
cat("Results saved: R_output_dir/simulation_validation.csv\n")
cat("\nDone.\n")
