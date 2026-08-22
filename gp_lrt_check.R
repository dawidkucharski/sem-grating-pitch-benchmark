# Reproduce GP LRT exactly with the same functions as gp_linearity.R
base_folder <- "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"
output_dir  <- file.path(base_folder, "R_output_dir")
mags <- c(400, 500, 750, 1000)

se_kernel <- function(X1, X2, sigma2, lengthscale) {
  sqdist <- outer(X1[,1], X2[,1], "-")^2 + outer(X1[,2], X2[,2], "-")^2
  sigma2 * exp(-0.5 * sqdist / lengthscale^2)
}
nlml <- function(par, X, y) {
  sigma2 <- exp(par[1]); lengthscale <- exp(par[2]); noise <- exp(par[3])
  n <- length(y)
  K <- se_kernel(X, X, sigma2, lengthscale) + noise * diag(n)
  L <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(L)) return(1e10)
  alpha <- backsolve(L, forwardsolve(L, y, upper.tri = FALSE))
  0.5 * t(y) %*% alpha + sum(log(diag(L))) + 0.5 * n * log(2 * pi)
}

for (mag in mags) {
  df <- read.csv(file.path(output_dir, sprintf("linearity_tiles_mag%d.csv", mag)))
  med <- median(df$local_pitch_um, na.rm = TRUE)
  madv <- mad(df$local_pitch_um, constant = 1.4826, na.rm = TRUE)
  keep <- abs(df$local_pitch_um - med) < 3 * madv
  df <- df[keep, ]
  n <- nrow(df)
  X <- cbind((df$tile_x - min(df$tile_x)) / max(1, max(df$tile_x) - min(df$tile_x)),
             (df$tile_y - min(df$tile_y)) / max(1, max(df$tile_y) - min(df$tile_y)))
  y <- df$local_pitch_um - mean(df$local_pitch_um)
  var0 <- mean(y^2)
  nll0 <- n/2 * (log(2*pi*var0) + 1)
  init <- c(log(var0), log(0.3), log(var0*0.1))
  opt <- optim(init, nlml, X = X, y = y, method = "L-BFGS-B",
               lower = c(log(1e-6), log(0.01), log(1e-8)),
               upper = c(log(1), log(2), log(1)),
               control = list(maxit = 500))
  lrt <- 2 * (nll0 - opt$value)
  pv <- pchisq(lrt, df = 2, lower.tail = FALSE)
  cat(sprintf("%d x: n=%d sigma2=%.3e ls=%.4f noise=%.3e nlml=%.3f nll0=%.3f LRT=%.2f p=%.3e\n",
      mag, n, exp(opt$par[1]), exp(opt$par[2]), exp(opt$par[3]), opt$value, nll0, lrt, pv))
}
