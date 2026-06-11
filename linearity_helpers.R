library(imager)
library(splines)
library(pracma)

base_folder <- "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"

config <- list(
  scale_per_pixel_um = 0.01333594,
  spline_spar        = 0.1
)

evaluate_linearity_from_tiles <- function(mag) {
  cat("\n==========\nLinearity from tiles for magnification:", mag, "\n==========\n")

  sed_folder <- file.path(base_folder, sprintf("mag%d_SED-BED", mag), "SED")
  sed_files <- list.files(sed_folder, full.names = TRUE,
                          pattern = "\\.(jpg|png|tif|bmp)$")
  if (length(sed_files) == 0) {
    cat("No SED images found for magnification", mag, "in", sed_folder, "\n")
    return(invisible(NULL))
  }

  tile_results <- list()
  for (f in sed_files) {
    img <- tryCatch({ load.image(f) }, error = function(e) NULL)
    if (is.null(img)) next

    if (dim(img)[3] > 1) {
      gray_img <- grayscale(img)
    } else {
      gray_img <- img
    }
    gray_mat <- as.matrix(gray_img[,,1])
    intensity_profile <- rowMeans(gray_mat)

    spline_fit <- smooth.spline(x = 1:length(intensity_profile),
                                y = intensity_profile,
                                spar = config$spline_spar)
    smoothed_profile <- spline_fit$y

    inverted_profile <- -smoothed_profile
    minima <- tryCatch({
      findpeaks(inverted_profile, nups = 1, ndowns = 1, threshold = 0.02)
    }, error = function(e) NULL)
    if (is.null(minima)) next
    minima_pos <- minima[,2]
    if (length(minima_pos) < 2) next

    dist_px <- diff(minima_pos)
    dist_um <- dist_px * config$scale_per_pixel_um
    local_pitch_um <- mean(dist_um, na.rm = TRUE)

    fname <- basename(f)
    m <- regexpr("X[0-9]+_Y[0-9]+", fname)
    tile_x <- NA_integer_
    tile_y <- NA_integer_
    if (m != -1) {
      xy_str <- regmatches(fname, m)
      parts <- strsplit(xy_str, "_")[[1]]
      x_str <- sub("X", "", parts[1])
      y_str <- sub("Y", "", parts[2])
      tile_x <- as.integer(x_str)
      tile_y <- as.integer(y_str)
    }

    tile_results[[length(tile_results) + 1]] <- data.frame(
      magnification = mag,
      image_path = f,
      tile_x = tile_x,
      tile_y = tile_y,
      local_pitch_um = local_pitch_um,
      stringsAsFactors = FALSE
    )
  }

  if (length(tile_results) == 0) {
    cat("No valid tile results for magnification", mag, "\n")
    return(invisible(NULL))
  }

  lin_df <- do.call(rbind, tile_results)

  has_x <- !is.na(lin_df$tile_x)
  has_y <- !is.na(lin_df$tile_y)

  if (any(has_x)) {
    lm_x <- lm(local_pitch_um ~ tile_x, data = lin_df[has_x,])
    cat("\nLinearity model along X (tiles):\n")
    print(summary(lm_x))
  } else {
    lm_x <- NULL
    cat("No valid tile_x indices; skipping X-linearity model.\n")
  }

  if (any(has_y)) {
    lm_y <- lm(local_pitch_um ~ tile_y, data = lin_df[has_y,])
    cat("\nLinearity model along Y (tiles):\n")
    print(summary(lm_y))
  } else {
    lm_y <- NULL
    cat("No valid tile_y indices; skipping Y-linearity model.\n")
  }

  output_dir <- file.path(base_folder, "R_output_dir")
  if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)
  lin_file <- file.path(output_dir, sprintf("linearity_tiles_mag%d.csv", mag))
  write.csv(lin_df, lin_file, row.names = FALSE)
  cat("Tile-based linearity data saved to:", lin_file, "\n")

  # Compute simple field non-linearity metrics (relative change across field)
  mean_pitch <- mean(lin_df$local_pitch_um, na.rm = TRUE)
  dx <- if (any(has_x)) max(lin_df$tile_x[has_x]) - min(lin_df$tile_x[has_x]) else NA_real_
  dy <- if (any(has_y)) max(lin_df$tile_y[has_y]) - min(lin_df$tile_y[has_y]) else NA_real_

  max_rel_change_x <- if (!is.null(lm_x) && !is.na(dx) && !is.na(mean_pitch) && mean_pitch != 0) {
    100 * abs(coef(lm_x)["tile_x"]) * dx / mean_pitch
  } else NA_real_

  max_rel_change_y <- if (!is.null(lm_y) && !is.na(dy) && !is.na(mean_pitch) && mean_pitch != 0) {
    100 * abs(coef(lm_y)["tile_y"]) * dy / mean_pitch
  } else NA_real_

  summary_row <- data.frame(
    magnification = mag,
    mean_pitch_um = mean_pitch,
    slope_x_um_per_index = if (!is.null(lm_x)) unname(coef(lm_x)["tile_x"]) else NA_real_,
    slope_y_um_per_index = if (!is.null(lm_y)) unname(coef(lm_y)["tile_y"]) else NA_real_,
    max_rel_change_x_pct = max_rel_change_x,
    max_rel_change_y_pct = max_rel_change_y,
    stringsAsFactors = FALSE
  )

  summary_file <- file.path(output_dir, "linearity_summary_overview.csv")
  if (file.exists(summary_file)) {
    old <- tryCatch(read.csv(summary_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(old)) {
      # replace row for this magnification if it exists
      old <- old[old$magnification != mag, ]
      summary_row <- rbind(old, summary_row)
    }
  }
  write.csv(summary_row, summary_file, row.names = FALSE)
  cat("Linearity overview updated in:", summary_file, "\n")

  plot_file <- file.path(base_folder, sprintf("linearity_tiles_mag%d.pdf", mag))
  pdf(plot_file, width = 7, height = 5)
  par(mfrow = c(1, 2))
  if (!is.null(lm_x) && any(has_x)) {
    plot(lin_df$tile_x[has_x], lin_df$local_pitch_um[has_x],
         xlab = "Tile Index X",
         ylab = "Local Groove Spacing (µm)",
         main = sprintf("Linearity Along X - mag %d", mag),
         pch = 16, col = rgb(0, 0, 1, 0.5))
    abline(lm_x, col = "red", lwd = 2)
  } else {
    plot.new(); title(main = "No X Data for Linearity Assessment")
  }
  if (!is.null(lm_y) && any(has_y)) {
    plot(lin_df$tile_y[has_y], lin_df$local_pitch_um[has_y],
         xlab = "Tile Index Y",
         ylab = "Local Groove Spacing (µm)",
         main = sprintf("Linearity Along Y - mag %d", mag),
         pch = 16, col = rgb(0, 0, 1, 0.5))
    abline(lm_y, col = "red", lwd = 2)
  } else {
    plot.new(); title(main = "No Y Data for Linearity Assessment")
  }
  dev.off()
  cat("Tile-based linearity plot saved to:", plot_file, "\n")

  invisible(lin_df)
}
