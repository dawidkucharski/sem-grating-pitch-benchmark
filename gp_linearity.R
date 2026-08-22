# =========================================================================
# gp_linearity.R
# Gaussian process regression for SEM field linearity assessment
# =========================================================================
# Replaces the independent X/Y OLS linear fits with a proper 2D spatial
# model. The squared-exponential kernel lengthscale quantifies the spatial
# scale of magnification variations: small lengthscale = localised
# distortions, large lengthscale = flat field (good linearity).
# =========================================================================

library(pbapply)

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------
base_folder <- "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"
output_dir  <- file.path(base_folder, "R_output_dir")
mags <- c(400, 500, 750, 1000)

# -------------------------------------------------------------------
# Squared-exponential (RBF) kernel
# -------------------------------------------------------------------
se_kernel <- function(X1, X2, sigma2, lengthscale) {
  # X1: n1 x 2 matrix, X2: n2 x 2 matrix
  # Returns: n1 x n2 covariance matrix
  sqdist <- outer(X1[,1], X2[,1], "-")^2 + outer(X1[,2], X2[,2], "-")^2
  sigma2 * exp(-0.5 * sqdist / lengthscale^2)
}

# -------------------------------------------------------------------
# Negative log marginal likelihood (to minimise)
# -------------------------------------------------------------------
nlml <- function(par, X, y) {
  # par = c(log_sigma2, log_lengthscale, log_noise)
  sigma2     <- exp(par[1])
  lengthscale <- exp(par[2])
  noise      <- exp(par[3])
  
  n <- length(y)
  K <- se_kernel(X, X, sigma2, lengthscale) + noise * diag(n)
  
  # Cholesky decomposition for stable computation
  L <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(L)) return(1e10)
  
  alpha <- backsolve(L, forwardsolve(L, y, upper.tri = FALSE))
  # log det K = 2 * sum(log(diag(L)))
  0.5 * t(y) %*% alpha + sum(log(diag(L))) + 0.5 * n * log(2 * pi)
}

# -------------------------------------------------------------------
# GP prediction at test points
# -------------------------------------------------------------------
gp_predict <- function(X_train, y_train, X_test, sigma2, lengthscale, noise) {
  n_train <- nrow(X_train)
  n_test  <- nrow(X_test)
  
  K_tt <- se_kernel(X_train, X_train, sigma2, lengthscale) + noise * diag(n_train)
  K_ts <- se_kernel(X_train, X_test,  sigma2, lengthscale)
  K_ss <- se_kernel(X_test,  X_test,  sigma2, lengthscale)
  
  L <- chol(K_tt)
  alpha <- backsolve(L, forwardsolve(L, y_train, upper.tri = FALSE))
  v <- forwardsolve(L, K_ts)
  
  mu_pred <- t(K_ts) %*% alpha
  # Predictive variance: diag(K_ss - t(v) %*% v)
  var_pred <- diag(K_ss) - colSums(v^2)
  var_pred[var_pred < 0] <- 0  # numerical safety
  
  list(mean = as.numeric(mu_pred), sd = sqrt(var_pred))
}

# -------------------------------------------------------------------
# Fit GP to tile-level linearity data for one magnification
# -------------------------------------------------------------------
fit_gp_linearity <- function(mag) {
  cat(sprintf("\n========== GP linearity: mag %d ==========\n", mag))
  
  # Load tile data
  csv_file <- file.path(output_dir, sprintf("linearity_tiles_mag%d.csv", mag))
  if (!file.exists(csv_file)) {
    cat(sprintf("  File not found: %s\n", csv_file))
    return(NULL)
  }
  
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d tiles\n", nrow(df)))
  
  # Remove outliers: pitch values more than 3 MAD from median
  med_pitch <- median(df$local_pitch_um, na.rm = TRUE)
  mad_pitch <- mad(df$local_pitch_um, constant = 1.4826, na.rm = TRUE)
  keep <- abs(df$local_pitch_um - med_pitch) < 3 * mad_pitch
  df <- df[keep, ]
  cat(sprintf("  After outlier removal: %d tiles\n", nrow(df)))
  
  # Normalise coordinates to [0,1] range for numerical stability
  tile_x_norm <- (df$tile_x - min(df$tile_x)) / max(1, max(df$tile_x) - min(df$tile_x))
  tile_y_norm <- (df$tile_y - min(df$tile_y)) / max(1, max(df$tile_y) - min(df$tile_y))
  
  X <- cbind(tile_x_norm, tile_y_norm)
  y <- df$local_pitch_um - mean(df$local_pitch_um)  # centre the response
  
  # Initial hyperparameters (on log scale)
  # sigma2: variance of pitch values
  # lengthscale: ~0.3 in normalised coordinates (moderate spatial correlation)
  # noise: ~10% of signal variance
  init_sigma2 <- var(y)
  init_par <- c(log(init_sigma2), log(0.3), log(init_sigma2 * 0.1))
  
  cat("  Optimising hyperparameters...\n")
  opt <- optim(init_par, nlml, X = X, y = y,
               method = "L-BFGS-B",
               lower = c(log(1e-6), log(0.01), log(1e-8)),
               upper = c(log(1), log(2), log(1)),
               control = list(maxit = 500))
  
  sigma2_opt     <- exp(opt$par[1])
  lengthscale_opt <- exp(opt$par[2])
  noise_opt      <- exp(opt$par[3])
  
  cat(sprintf("  sigma^2 = %.6f\n", sigma2_opt))
  cat(sprintf("  lengthscale = %.4f (normalised coord)\n", lengthscale_opt))
  cat(sprintf("  noise variance = %.8f\n", noise_opt))
  cat(sprintf("  NLML = %.3f\n", opt$value))
  
  # Interpret lengthscale in physical units
  # In normalised coords, lengthscale=0.3 means correlation decays over ~30% of the field
  # Convert to physical units using the stage-coordinate field extent
  # (Table: linearity); these are the true physical positions of the tile centres.
  stage_field_width_um <- c("400" = 305, "500" = 237, "750" = 153, "1000" = 119)[as.character(mag)]
  lengthscale_um <- lengthscale_opt * stage_field_width_um
  
  cat(sprintf("  Approx. physical lengthscale: %.1f µm\n", lengthscale_um))
  cat(sprintf("  Stage field width: %.0f µm\n", stage_field_width_um))
  
  # Predict on a fine grid for visualisation
  ngrid <- 50
  grid_x <- seq(0, 1, length.out = ngrid)
  grid_y <- seq(0, 1, length.out = ngrid)
  X_grid <- as.matrix(expand.grid(grid_x, grid_y))
  
  pred <- gp_predict(X, y, X_grid, sigma2_opt, lengthscale_opt, noise_opt)
  
  # Denormalise predictions
  pred_mean <- pred$mean + mean(df$local_pitch_um)
  pred_sd   <- pred$sd
  
  # Maximum relative distortion across the field (normalised)
  max_distortion_pct <- 100 * (max(pred$mean) - min(pred$mean)) / mean(df$local_pitch_um)
  cat(sprintf("  Max GP distortion range: %.3f%% of mean pitch\n", max_distortion_pct))
  
  # ---- Generate figure ----
  pdf(file.path(base_folder, sprintf("gp_linearity_mag%d.pdf", mag)),
      width = 10, height = 5)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # Distortion map
  z_mean <- matrix(pred_mean, nrow = ngrid, ncol = ngrid)
  image(grid_x, grid_y, z_mean,
        col = hcl.colors(64, "RdBu", rev = TRUE),
        xlab = "Normalised X", ylab = "Normalised Y",
        main = sprintf("GP Distortion Map — mag %d×", mag))
  contour(grid_x, grid_y, z_mean, add = TRUE, col = "grey30", lwd = 0.5)
  points(tile_x_norm, tile_y_norm, pch = ".", col = "grey20")
  
  # Uncertainty map (predictive SD)
  z_sd <- matrix(pred_sd, nrow = ngrid, ncol = ngrid)
  image(grid_x, grid_y, z_sd,
        col = hcl.colors(64, "YlOrRd", rev = TRUE),
        xlab = "Normalised X", ylab = "Normalised Y",
        main = sprintf("GP Predictive SD — mag %d×", mag))
  points(tile_x_norm, tile_y_norm, pch = ".", col = "grey80")
  
  dev.off()
  cat(sprintf("  Figure saved: gp_linearity_mag%d.pdf\n", mag))
  
  # Return summary
  data.frame(
    magnification      = mag,
    n_tiles            = nrow(df),
    sigma2             = sigma2_opt,
    lengthscale_norm   = lengthscale_opt,
    lengthscale_um     = lengthscale_um,
    noise_variance     = noise_opt,
    max_distortion_pct = max_distortion_pct,
    nlml               = opt$value,
    stringsAsFactors   = FALSE
  )
}

# =========================================================================
# Run for all magnifications
# =========================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("GP Linearity Analysis\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

gp_summaries <- list()
for (mag in mags) {
  result <- fit_gp_linearity(mag)
  if (!is.null(result)) {
    gp_summaries[[length(gp_summaries) + 1]] <- result
  }
}

gp_summary_df <- do.call(rbind, gp_summaries)

# Save summary
summary_file <- file.path(output_dir, "gp_linearity_summary.csv")
write.csv(gp_summary_df, summary_file, row.names = FALSE)
cat(sprintf("\nSummary saved: %s\n", summary_file))

# Print overview
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("GP Linearity Summary\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
print(gp_summary_df[, c("magnification", "n_tiles", "lengthscale_um",
                          "max_distortion_pct")],
      row.names = FALSE)

# ---- Combined figure: all four mags ----
pdf(file.path(base_folder, "gp_linearity_combined.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (mag in mags) {
  csv_file <- file.path(output_dir, sprintf("linearity_tiles_mag%d.csv", mag))
  if (!file.exists(csv_file)) next
  
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  med_pitch <- median(df$local_pitch_um, na.rm = TRUE)
  mad_pitch <- mad(df$local_pitch_um, constant = 1.4826, na.rm = TRUE)
  df <- df[abs(df$local_pitch_um - med_pitch) < 3 * mad_pitch, ]
  
  tile_x_norm <- (df$tile_x - min(df$tile_x)) / max(1, max(df$tile_x) - min(df$tile_x))
  tile_y_norm <- (df$tile_y - min(df$tile_y)) / max(1, max(df$tile_y) - min(df$tile_y))
  
  X <- cbind(tile_x_norm, tile_y_norm)
  y <- df$local_pitch_um - mean(df$local_pitch_um)
  
  init_par <- c(log(var(y)), log(0.3), log(var(y) * 0.1))
  opt <- optim(init_par, nlml, X = X, y = y,
               method = "L-BFGS-B",
               lower = c(log(1e-6), log(0.01), log(1e-8)),
               upper = c(log(1), log(2), log(1)),
               control = list(maxit = 500))
  
  sigma2_opt <- exp(opt$par[1]); ls_opt <- exp(opt$par[2]); n_opt <- exp(opt$par[3])
  
  ngrid <- 40
  grid_x <- seq(0, 1, length.out = ngrid)
  grid_y <- seq(0, 1, length.out = ngrid)
  X_grid <- as.matrix(expand.grid(grid_x, grid_y))
  
  pred <- gp_predict(X, y, X_grid, sigma2_opt, ls_opt, n_opt)
  z <- matrix(pred$mean + mean(df$local_pitch_um), nrow = ngrid, ncol = ngrid)
  
  # Use a consistent colour range centred on the median pitch for this mag
  pitch_range <- range(z, na.rm = TRUE)
  stage_field_width_um <- c("400" = 305, "500" = 237, "750" = 153, "1000" = 119)[as.character(mag)]
  image(grid_x, grid_y, z,
        col = hcl.colors(64, "RdBu", rev = TRUE),
        zlim = pitch_range,
        xlab = "Normalised X", ylab = "Normalised Y",
        main = sprintf("%d×  (lengthscale = %.0f µm)", mag,
                       ls_opt * stage_field_width_um))
  contour(grid_x, grid_y, z, add = TRUE, col = "grey30", lwd = 0.5)
}

dev.off()
cat(sprintf("Combined figure: gp_linearity_combined.pdf\n"))

cat("\nDone.\n")
