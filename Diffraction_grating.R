# Load necessary libraries
library(imager)    # For image processing
library(splines)   # For B-spline smoothing
library(pracma)    # For peak detection
library(parallel)  # For parallel processing
library(lubridate) # For timing
library(progress)  # For progress bar
library(pbapply)   # For progress bar in parallel

# -------------------------------------------------------------------
# Global configuration (algorithm parameters, paths, constants)
# -------------------------------------------------------------------

base_folder <- "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"

config <- list(
  scale_per_pixel_um      = 0.01333594,   # pixel size [um]
  min_groove_fraction     = 0.6,          # reject if any spacing < fraction * mean spacing
  spline_spar             = 0.1,          # smoothing parameter for smooth.spline
  roughness_decimals      = 3,            # rounding for Ra/Rz etc.
  k_coverage              = 3,            # coverage factor for expanded uncertainty
  typeB_uncertainty_um    = 0.01333594,   # type B component [um]
  angle_window_deg        = 3,            # window around dominant angle for core lines
  grad_quantile_mag500    = 0.70,         # gradient threshold quantile for mag=500
  grad_quantile_other     = 0.80          # gradient threshold quantile for other mags
)


# -------------------------------------------------------------------
# Optional: Linearity evaluation on stitched images
# -------------------------------------------------------------------

# This helper evaluates local pitch (groove spacing) as a function of
# position along the stitching direction, using one stitched image per
# magnification. It does not modify the main analysis; it only generates
# additional CSV and an optional plot.

evaluate_linearity_for_mag <- function(mag,
                                       stitched_path = NULL,
                                       n_profiles = 20,
                                       profile_height_px = 20) {
  cat("\n==========\nLinearity evaluation for magnification:", mag, "\n==========\n")

  # If no stitched image path is supplied, try a default convention
  if (is.null(stitched_path)) {
    # Example: mag1000/stitched/mag1000_stitched.png
    stitched_path <- file.path(base_folder,
                               sprintf("mag%d", mag),
                               "stitched",
                               sprintf("mag%d_stitched.png", mag))
  }

  if (!file.exists(stitched_path)) {
    cat("Stitched image not found for magnification", mag, "at:\n",
        stitched_path, "\n")
    return(invisible(NULL))
  }

  cat("Using stitched image:\n", stitched_path, "\n")

  img <- tryCatch({
    load.image(stitched_path)
  }, error = function(e) {
    cat("Error loading stitched image:", e$message, "\n")
    return(NULL)
  })

  if (is.null(img)) return(invisible(NULL))

  # Convert to grayscale
  if (dim(img)[3] > 1) {
    gray_img <- grayscale(img)
  } else {
    gray_img <- img
  }

  # Work in matrix form
  gray_mat <- as.matrix(gray_img[,,1])
  n_rows <- nrow(gray_mat)
  n_cols <- ncol(gray_mat)

  # Choose profile center rows uniformly across the image height
  if (profile_height_px > n_rows) profile_height_px <- n_rows
  row_centers <- round(seq(from = floor(profile_height_px/2) + 1,
                           to   = n_rows - floor(profile_height_px/2),
                           length.out = n_profiles))

  results <- list()

  for (rc in row_centers) {
    row_start <- max(1, rc - floor(profile_height_px/2))
    row_end   <- min(n_rows, rc + floor(profile_height_px/2))

    sub_mat <- gray_mat[row_start:row_end, , drop = FALSE]
    # intensity profile along columns (mean over selected rows)
    intensity_profile <- colMeans(sub_mat)

    # smooth profile with the same spline approach
    spline_fit <- smooth.spline(x = 1:length(intensity_profile),
                                y = intensity_profile,
                                spar = config$spline_spar)
    smoothed_profile <- spline_fit$y

    # detect minima via inverted profile
    inverted_profile <- -smoothed_profile
    minima <- tryCatch({
      findpeaks(inverted_profile, nups = 1, ndowns = 1, threshold = 0.02)
    }, error = function(e) {
      NULL
    })

    if (is.null(minima)) {
      next
    }

    minima_pos <- minima[, 2]
    if (length(minima_pos) < 2) next

    # distances between successive minima in pixels and um
    dist_px <- diff(minima_pos)
    dist_um <- dist_px * config$scale_per_pixel_um

    local_pitch_um <- mean(dist_um, na.rm = TRUE)
    center_x_px <- mean(range(minima_pos))
    center_x_um <- center_x_px * config$scale_per_pixel_um

    results[[length(results) + 1]] <- data.frame(
      magnification = mag,
      row_center_px = rc,
      x_pos_px = center_x_px,
      x_pos_um = center_x_um,
      local_pitch_um = local_pitch_um
    )
  }

  if (length(results) == 0) {
    cat("No valid profiles for linearity evaluation at magnification", mag, "\n")
    return(invisible(NULL))
  }

  lin_df <- do.call(rbind, results)

  # Fit a simple linear model: pitch vs position
  lm_fit <- lm(local_pitch_um ~ x_pos_um, data = lin_df)
  cat("\nLinearity model (local pitch vs position):\n")
  print(summary(lm_fit))

  # Save CSV
  output_dir <- file.path(base_folder, "R_output_dir")
  if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)

  lin_file <- file.path(output_dir, sprintf("linearity_mag%d.csv", mag))
  write.csv(lin_df, lin_file, row.names = FALSE)
  cat("Linearity data saved to:", lin_file, "\n")

  # Optional quick plot to PDF
  plot_file <- file.path(base_folder, sprintf("linearity_mag%d.pdf", mag))
  pdf(plot_file, width = 7, height = 5)
  plot(lin_df$x_pos_um, lin_df$local_pitch_um,
      xlab = "Position Across Field (µm)",
      ylab = "Local Groove Spacing (µm)",
      main = sprintf("Magnification Linearity - mag %d", mag),
       pch = 16, col = rgb(0, 0, 1, 0.5))
  abline(lm_fit, col = "red", lwd = 2)
  dev.off()
  cat("Linearity plot saved to:", plot_file, "\n")

  invisible(lin_df)
}

run_analysis_for_mag <- function(mag) {
  cat("\n==========\nAnalysis for magnification:", mag, "\n==========\n")

  sed_folder <- file.path(base_folder, sprintf("mag%d_SED-BED", mag), "SED")
  beds_folder <- file.path(base_folder, sprintf("mag%d_SED-BED", mag), "BED-S")

  # Get list of images from both folders
  sed_files <- list.files(sed_folder, full.names = TRUE, pattern = "\\.(jpg|png|tif|bmp)$")
  beds_files <- list.files(beds_folder, full.names = TRUE, pattern = "\\.(jpg|png|tif|bmp)$")

  # Ensure correct pairing by names
  sed_names <- sub("_SED_", "", tools::file_path_sans_ext(basename(sed_files)))
  beds_names <- sub("_BED-S_", "", tools::file_path_sans_ext(basename(beds_files)))
  common_names <- intersect(sed_names, beds_names)

  if (length(common_names) == 0) {
    stop("No matching image names found between SED and BED-S folders for magnification ", mag)
  }

  # Filter only paired files
  sed_files <- sed_files[sed_names %in% common_names]
  beds_files <- beds_files[beds_names %in% common_names]

  # Randomly sample the required proportions
  set.seed(123)  # For reproducibility
  n_beds <- floor(0.5 * length(common_names))
  n_sed <- length(common_names) - n_beds

  if (n_beds > length(beds_files) || n_sed > length(sed_files)) {
    stop("Sample size exceeds available images for magnification ", mag)
  }

  selected_beds <- sample(beds_files, n_beds)
  selected_sed <- sample(sed_files, n_sed)
  mixed_files <- c(selected_beds, selected_sed)


  # helper to tag detector type based on file name
  detect_detector <- function(path) {
    name <- basename(path)
    if (grepl("_SED_", name, ignore.case = TRUE)) return("SED")
    if (grepl("_BED-S_", name, ignore.case = TRUE)) return("BED-S")
    return("UNKNOWN")
  }

  # container for rejected images
  rejected <- list()

  process_image <- function(image_path, mean_intensity, pb) {
  # Load the image with error handling
  image <- tryCatch({
    load.image(image_path)
  }, error = function(e) {
    cat("Error loading image:", image_path, "\n")
    return(NULL)
  })
  
  if (is.null(image)) return(NULL)
  
  # Ensure the image is grayscale
  if (dim(image)[3] > 1) {
    gray_image <- grayscale(image)
  } else {
    gray_image <- image
  }
  
  # Normalize image intensity to match the overall mean intensity
  current_mean <- mean(gray_image)
  gray_image <- gray_image + (mean_intensity - current_mean)
  
  # Convert the grayscale image to a 2D matrix
  gray_matrix <- as.matrix(gray_image[,,1])  # Ensure only 2D data
  
  # Aggregate pixel intensities along rows
  intensity_profile <- rowMeans(gray_matrix)
  
  # Smooth the intensity profile using B-splines
  spline_fit <- smooth.spline(x = 1:length(intensity_profile), y = intensity_profile,
                              spar = config$spline_spar)
  smoothed_profile <- spline_fit$y
  
  # Detect minima in the smoothed profile (invert the profile to find minima as peaks)
  inverted_profile <- -smoothed_profile
  minima <- tryCatch({
    findpeaks(inverted_profile, nups = 1, ndowns = 1, threshold = 0.02)
  }, error = function(e) {
    return(NULL)
  })
  
  # Detect maxima directly from the smoothed profile
  maxima <- tryCatch({
    findpeaks(smoothed_profile, nups = 1, ndowns = 1, threshold = 0.02)
  }, error = function(e) {
    return(NULL)
  })
  
  # Return NULL if no minima or maxima were detected
  if (is.null(minima) || is.null(maxima)) {
    return(NULL)  
  }
  
  # Extract minima and maxima positions
  minima_positions <- minima[, 2]
  maxima_positions <- maxima[, 2]
  
  # Calculate distances between minima
  distances <- diff(minima_positions)
  
  # Convert distances to real-world units (micrometers)
  real_distances <- distances * config$scale_per_pixel_um
  
  # Calculate mean distance for the current image
  mean_distance <- round(mean(real_distances, na.rm = TRUE), 3)
  
  # Filter out results where any groove distance is smaller than fraction of the average distance
  if (any(real_distances < (mean_distance * config$min_groove_fraction))) {
    rejected[[length(rejected) + 1]] <<- list(
      image_path = image_path,
      reason = "too_small_spacing",
      mean_distance_um = mean_distance
    )
    return(NULL)  # Exclude this result if the condition is met
  }
  
  # Calculate the mean height between minima and maxima
  mean_heights <- sapply(1:length(minima_positions), function(i) {
    if (i <= length(maxima_positions)) {
      min_pos <- minima_positions[i]
      max_pos <- maxima_positions[i]
      if (max_pos > min_pos) {
        return(mean(smoothed_profile[min_pos:max_pos]))
      }
    }
    return(NA)
  })
  
  # Calculate the mean height for this image
  mean_height <- round(mean(mean_heights, na.rm = TRUE), 3)

  # Intensity-profile "roughness" (not true height Ra/Rz):
  # remove linear trend from the smoothed profile and compute
  # Ra_int = mean(|residuals|), Rz_int = max(residuals) - min(residuals)
  x_idx <- 1:length(smoothed_profile)
  profile_lm <- lm(smoothed_profile ~ x_idx)
  profile_resid <- residuals(profile_lm)
  Ra <- mean(abs(profile_resid), na.rm = TRUE)
  profile_max <- max(profile_resid, na.rm = TRUE)
  profile_min <- min(profile_resid, na.rm = TRUE)
  Rz <- profile_max - profile_min
  
  # Update progress bar
  if (!is.null(pb)) {
    pb$tick()
  }
  
  return(list(mean_distance = mean_distance, mean_height = mean_height,
              real_distances = real_distances, image_path = image_path,
              detector = detect_detector(image_path),
              intensity_profile = intensity_profile, smoothed_profile = smoothed_profile,
              minima_positions = minima_positions, maxima_positions = maxima_positions,
              Ra = Ra, Rz = Rz))

  }

  # Step 3: Parallel processing with progress bar using pblapply
  library(pbapply)  # Make sure pbapply is loaded

  num_cores <- detectCores() - 1
  start_time <- Sys.time()

  # Create the progress bar
  pb <- progress_bar$new(
    format = "  [:bar] :percent :elapsed",
    total = length(mixed_files),
    clear = FALSE, width = 60
  )

  # Calculate overall mean intensity in parallel
  mean_intensity_all <- mean(unlist(pblapply(mixed_files, function(image_path) {
    image <- tryCatch({
      load.image(image_path)
    }, error = function(e) {
      return(NULL)
    })
    if (!is.null(image)) {
      if (dim(image)[3] > 1) {
        gray_image <- grayscale(image)
      } else {
        gray_image <- image
      }
      return(mean(gray_image))
    } else {
      return(NA)
    }
  }, cl = num_cores)))  # Specify number of cores to use

  # Check for any NA values (in case some images failed to load)
  if (any(is.na(mean_intensity_all))) {
    cat("Some images failed to load or process. These will be excluded from further analysis.\n")
    mean_intensity_all <- mean(mean_intensity_all, na.rm = TRUE)
  }

  # Process images in parallel using pblapply with progress bar
  results <- pblapply(mixed_files, process_image, mean_intensity = mean_intensity_all, pb = pb, cl = num_cores)
  results <- results[!sapply(results, is.null)]

  # Step 4: Calculate and print results
  if (length(results) > 0) {
  all_mean_distances <- sapply(results, function(x) x$mean_distance)
  all_mean_heights <- sapply(results, function(x) x$mean_height)
  all_Ra <- sapply(results, function(x) x$Ra)
  all_Rz <- sapply(results, function(x) x$Rz)
  
  # Normality Test for mean distances
  shapiro_test <- shapiro.test(all_mean_distances)
  is_normal <- shapiro_test$p.value > 0.05
  
  # Overall Metrics based on distribution type
  if (is_normal) {
    overall_mean_distance <- round(mean(all_mean_distances, na.rm = TRUE), 3)
    overall_mean_height <- round(mean(all_mean_heights, na.rm = TRUE), 3)
    type_a_uncertainty <- round(sd(all_mean_distances, na.rm = TRUE) / sqrt(length(all_mean_distances)), 3)
  } else {
    overall_mean_distance <- round(median(all_mean_distances, na.rm = TRUE), 3)
    overall_mean_height <- round(median(all_mean_heights, na.rm = TRUE), 3)
    mad_mean_distances <- mad(all_mean_distances, constant = 1.4826, na.rm = TRUE)  # Robust estimator of variability
    type_a_uncertainty <- round(mad_mean_distances / sqrt(length(all_mean_distances)), 3)
  }
  
  # Type B Uncertainty (provided)
  type_b_uncertainty <- config$typeB_uncertainty_um
  
  # Combined and Expanded Uncertainty
  combined_uncertainty <- round(sqrt(type_a_uncertainty^2 + type_b_uncertainty^2), 3)
  expanded_uncertainty <- round(config$k_coverage * combined_uncertainty, 3)
  
  # Overall Roughness Parameters (mean calculation, normality assumption isn't critical)
  overall_Ra <- round(mean(all_Ra, na.rm = TRUE), config$roughness_decimals)
  overall_Rz <- round(mean(all_Rz, na.rm = TRUE), config$roughness_decimals)
  
  # Output Results
  cat("\nDistribution Type for Groove Spacing: ", ifelse(is_normal, "Normal", "Non-Normal"), "\n")
  cat("\nOverall Groove Spacing (µm):", sprintf("%.3f ± %.3f", overall_mean_distance, expanded_uncertainty), "\n")
  cat("\nOverall Height (between minima and maxima):", sprintf("%.3f ± %.3f", overall_mean_height, expanded_uncertainty), "\n")  
  cat("Roughness Parameter Ra (µm):", sprintf("%.3f", overall_Ra), "\n")
  cat("Roughness Parameter Rz (µm):", sprintf("%.3f", overall_Rz), "\n")
  
  # Magnification Error
  delta_M <- ((overall_mean_distance - (1000 / 600)) / (1000 / 600)) * 100  
  cat("\nMagnification error:", sprintf("%.3f%%", abs(delta_M)), "\n")
  
  # Error in Groove Spacing
  Real_d <- round((1 / 600 * 1000), 3)
  Error <- abs(overall_mean_distance - Real_d)
  cat("\nError in Groove Spacing (µm):", sprintf("%.3f", Error), "\n")


  end_time <- Sys.time()
  
  # Time calculation
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  hours <- as.integer(time_taken %/% 3600)
  minutes <- as.integer((time_taken %% 3600) %/% 60)
  seconds <- as.integer(time_taken %% 60)
  
  cat("\nTime Taken (h:m:s):", sprintf("%02d:%02d:%02d", hours, minutes, seconds), "\n")
  
  # Count number of images evaluated
  num_images <- length(results)
  cat("\nNumber of images evaluated:", num_images, "\n")

  # Count number of rejected images
  n_rejected <- length(rejected)
  if (n_rejected > 0) {
    cat("Number of rejected images:", n_rejected, "\n")
  }

  # Export summary metrics to CSV in R_output
  # Use a directory named R_output_dir to avoid clashing with
  # an existing file named 'R_output' in the base folder.
  output_dir <- file.path(base_folder, "R_output_dir")
  if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)

  # split summary also by detector (SED vs BED-S)
  det_vec <- sapply(results, function(x) x$detector)
  sed_idx <- which(det_vec == "SED")
  beds_idx <- which(det_vec == "BED-S")

  summary_by_det <- function(idx) {
    if (length(idx) == 0) return(data.frame())
    d <- sapply(results[idx], function(x) x$mean_distance)
    h <- sapply(results[idx], function(x) x$mean_height)
    Ra_det <- sapply(results[idx], function(x) x$Ra)
    Rz_det <- sapply(results[idx], function(x) x$Rz)
    if (is_normal) {
      d_mean <- round(mean(d, na.rm = TRUE), 3)
      h_mean <- round(mean(h, na.rm = TRUE), 3)
      type_a_det <- round(sd(d, na.rm = TRUE) / sqrt(length(d)), 3)
    } else {
      d_mean <- round(median(d, na.rm = TRUE), 3)
      h_mean <- round(median(h, na.rm = TRUE), 3)
      mad_d <- mad(d, constant = 1.4826, na.rm = TRUE)
      type_a_det <- round(mad_d / sqrt(length(d)), 3)
    }
    combined_det <- round(sqrt(type_a_det^2 + type_b_uncertainty^2), 3)
    expanded_det <- round(config$k_coverage * combined_det, 3)
    Ra_mean_det <- round(mean(Ra_det, na.rm = TRUE), config$roughness_decimals)
    Rz_mean_det <- round(mean(Rz_det, na.rm = TRUE), config$roughness_decimals)
    delta_M_det <- ((d_mean - (1000 / 600)) / (1000 / 600)) * 100
    Error_det <- abs(d_mean - Real_d)
    data.frame(
      magnification = mag,
      detector = if (length(idx) > 0 && det_vec[idx[1]] == "SED") "SED" else "BED-S",
      n_images = length(idx),
      groove_spacing_um = d_mean,
      groove_spacing_unc_um = expanded_det,
      height_um = h_mean,
      height_unc_um = expanded_det,
      Ra_um = Ra_mean_det,
      Rz_um = Rz_mean_det,
      magnification_error_pct = abs(delta_M_det),
      groove_spacing_error_um = Error_det,
      stringsAsFactors = FALSE
    )
  }

  groove_summary_overall <- data.frame(
    magnification = mag,
    n_images = num_images,
    groove_spacing_um = overall_mean_distance,
    groove_spacing_unc_um = expanded_uncertainty,
    height_um = overall_mean_height,
    height_unc_um = expanded_uncertainty,
    Ra_um = overall_Ra,
    Rz_um = overall_Rz,
    magnification_error_pct = abs(delta_M),
    groove_spacing_error_um = Error,
    stringsAsFactors = FALSE
  )

  groove_summary_sed  <- summary_by_det(sed_idx)
  groove_summary_beds <- summary_by_det(beds_idx)

  groove_file <- file.path(output_dir, sprintf("groove_summary_mag%d.csv", mag))
  write.csv(groove_summary_overall, groove_file, row.names = FALSE)
  cat("\nGroove summary saved to:", groove_file, "\n")

  if (nrow(groove_summary_sed) > 0) {
    groove_file_sed <- file.path(output_dir, sprintf("groove_summary_mag%d_SED.csv", mag))
    write.csv(groove_summary_sed, groove_file_sed, row.names = FALSE)
    cat("SED groove summary saved to:", groove_file_sed, "\n")
  }
  if (nrow(groove_summary_beds) > 0) {
    groove_file_beds <- file.path(output_dir, sprintf("groove_summary_mag%d_BED-S.csv", mag))
    write.csv(groove_summary_beds, groove_file_beds, row.names = FALSE)
    cat("BED-S groove summary saved to:", groove_file_beds, "\n")
  }

  # export rejected images, if any
  if (length(rejected) > 0) {
    rej_df <- do.call(rbind, lapply(rejected, as.data.frame))
    rej_file <- file.path(output_dir, sprintf("rejected_images_mag%d.csv", mag))
    write.csv(rej_df, rej_file, row.names = FALSE)
    cat("Rejected images saved to:", rej_file, "\n")
  }
  # Step 4: Plotting results with original, SED, and BED-S images
  if (length(results) > 0 && interactive()) {
    # Get the first result for the profile and images
    first_result <- results[[1]]
    intensity_profile <- first_result$intensity_profile
    smoothed_profile <- first_result$smoothed_profile
    minima_positions <- first_result$minima_positions
    maxima_positions <- first_result$maxima_positions
    # Try to load the original image; if it fails, skip plotting
    original_image <- tryCatch({
      load.image(first_result$image_path)
    }, error = function(e) {
      cat("Could not load original image for plotting:\n", first_result$image_path, "\n")
      return(NULL)
    })

    if (is.null(original_image)) {
      cat("Skipping image-based plots due to missing original image.\n")
    } else {
      # Extract the name of the current image (for matching SED and BED-S)
      current_image_name <- tools::file_path_sans_ext(basename(first_result$image_path))
    
      # Remove common prefixes or suffixes for better matching (e.g., '_SED_' and '_BED-S_')
      core_name <- sub("_SED_|_BED-S_", "", current_image_name)
    
      # Find the corresponding SED and BED-S images by name
      sed_image_path <- sed_files[grep(core_name, sed_names)]
      bed_image_path <- beds_files[grep(core_name, beds_names)]
    
      # Check if corresponding images exist
      if (length(sed_image_path) == 0 | length(bed_image_path) == 0) {
        cat("No matching SED or BED-S images found for", core_name, "\n")
      } else {
        # Load the corresponding SED and BED-S images with error handling
        sed_image <- tryCatch({
          load.image(sed_image_path)
        }, error = function(e) {
          cat("Could not load SED image for plotting:\n", sed_image_path, "\n")
          return(NULL)
        })

        bed_image <- tryCatch({
          load.image(bed_image_path)
        }, error = function(e) {
          cat("Could not load BED-S image for plotting:\n", bed_image_path, "\n")
          return(NULL)
        })

        if (is.null(sed_image) || is.null(bed_image)) {
          cat("Skipping SED/BED-S image plots due to missing images.\n")
        } else {
          # Set up the plotting area (2x2 grid)
          par(mfrow = c(2, 2), mar = c(5, 5, 2, 1))
          
          # Plot the intensity profile with smoothing
          plot(1:length(intensity_profile), intensity_profile, type = "l", col = "gray",
               main = "Intensity Profile with B-Spline Smoothing",
               xlab = "Row Index", ylab = "Mean Intensity", lwd = 2)
          lines(1:length(smoothed_profile), smoothed_profile, col = "blue", lwd = 2)
          points(minima_positions, smoothed_profile[minima_positions], col = "red", pch = 19)
          points(maxima_positions, smoothed_profile[maxima_positions], col = "green", pch = 19)
          
          # Plot the original combined image
          plot(original_image, main = paste("Original Combined Image: ", core_name))
          
          # Plot the corresponding SED image
          plot(sed_image, main = paste("Corresponding SED Image: ", core_name))
          
          # Plot the corresponding BED-S image
          plot(bed_image, main = paste("Corresponding BED-S Image: ", core_name))
        }
      }
    }
  } else {
    cat("\nNo valid results found.\n")
  }
}


# -------------------------------------------------------------------
# Linearity evaluation based on per-image mean spacing across tiles
# -------------------------------------------------------------------

# This helper uses per-image mean groove spacing for all SED tiles
# of a given magnification. It parses tile indices from filenames
# (e.g. ..._SED_X003_Y014.png) and fits simple models of local pitch
# versus tile X/Y index to assess linearity across the stitched area.

evaluate_linearity_from_tiles <- function(mag) {
  cat("\n==========\nLinearity from tiles for magnification:", mag, "\n==========\n")

  sed_folder <- file.path(base_folder, sprintf("mag%d_SED-BED", mag), "SED")
  sed_files <- list.files(sed_folder, full.names = TRUE, pattern = "\\.(jpg|png|tif|bmp)$")
  if (length(sed_files) == 0) {
    cat("No SED images found for magnification", mag, "in", sed_folder, "\n")
    return(invisible(NULL))
  }

  # For each SED image, reuse the same processing logic in a simplified form
  tile_results <- list()
  for (f in sed_files) {
    img <- tryCatch({ load.image(f) }, error = function(e) NULL)
    if (is.null(img)) next

    # grayscale
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

    # parse tile indices Xnnn_Ymmm from filename
    fname <- basename(f)
    m <- regexpr("X[0-9]+_Y[0-9]+", fname)
    tile_x <- NA_integer_
    tile_y <- NA_integer_
    if (m != -1) {
      xy_str <- regmatches(fname, m)
      # xy_str like "X003_Y014"
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

  # Fit linear models if we have non-NA indices
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

  # Save CSV
  output_dir <- file.path(base_folder, "R_output_dir")
  if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)
  lin_file <- file.path(output_dir, sprintf("linearity_tiles_mag%d.csv", mag))
  write.csv(lin_df, lin_file, row.names = FALSE)
  cat("Tile-based linearity data saved to:", lin_file, "\n")

  # Optional plots
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


#-------------------Parallelism_evaluation------------------------
# This block evaluates parallelism of grooves for SED images in the
# current magnification folder. For robustness (especially at mag=500),
# it aggregates lines from all SED images instead of a single one.

if (length(sed_files) > 0) {
  all_fit_lines <- list()
  for (pe_image_path in sed_files) {
    cat("\n[Parallelism evaluation] Using image:\n", pe_image_path, "\n")
    img <- load.image(pe_image_path)
    channels <- spectrum(img)
    cat("Number of channels:", channels, "\n")
    if (channels == 1) {
      gray_img <- img
    } else if (channels >= 3) {
      R_ch <- R(img)
      G_ch <- G(img)
      B_ch <- B(img)
      gray_img <- 0.299 * R_ch + 0.587 * G_ch + 0.114 * B_ch
    } else {
      stop("Unsupported number of channels")
    }
    # 1. Gaussian smoothing
    blurred <- isoblur(gray_img, sigma = 9.0)
    # 2. Sobel gradients
    gx <- imgradient(blurred, "x")
    gy <- imgradient(blurred, "y")
    grad_mag <- sqrt(gx^2 + gy^2)
    # 3. Thresholds based on configuration
    if (mag == 500) {
      threshold <- quantile(grad_mag, config$grad_quantile_mag500)
    } else {
      threshold <- quantile(grad_mag, config$grad_quantile_other)
    }
    edges <- grad_mag > threshold
    # 3. Segment connected edge regions (label)
    labelled <- label(edges)
    # 4. Convert to data frame and filter edge pixels
    df <- as.data.frame(labelled) %>% dplyr::filter(value > 0)
    # 5. Fit a line to each component (group = value)
    fit_lines <- df %>%
      dplyr::group_by(value) %>%
      # Lower the pixel-count threshold to approach ~9 lines
      dplyr::filter(dplyr::n() > if (mag == 500) 80 else 200) %>%
      dplyr::summarise(
        fit = list(lm(y ~ x)),
        angle_deg = atan(coef(lm(y ~ x))[2]) * 180 / pi,
        .groups = "drop"
      ) %>%
      # Allow a small deviation from 90°
      dplyr::filter(abs((angle_deg %% 180) - 90) < 5)

    # Collect fitted lines for this image
    all_fit_lines[[pe_image_path]] <- fit_lines
  }

  # Bind all fitted lines from all SED images
  fit_lines_all <- dplyr::bind_rows(all_fit_lines)

  # 6. Determine dominant groove direction and robust parallelism error
  # 6a. Dominant angle (median is robust against outliers)
  dominant_angle <- stats::median(fit_lines_all$angle_deg, na.rm = TRUE)

  # 6b. Keep only lines close to dominant direction (± angle_window_deg)
  fit_lines_core <- fit_lines_all %>%
    dplyr::filter(abs(angle_deg - dominant_angle) < config$angle_window_deg)

  # 6c. Compute deviations from dominant angle
  angle_dev <- fit_lines_core$angle_deg - dominant_angle

  # 6d. Robust sigma using MAD of deviations
  mad_dev <- stats::mad(angle_dev, constant = 1.4826, na.rm = TRUE)
  sigma_parallel <- mad_dev

  # 6e. Deviation from 90 degrees (global tilt)
  delta_alpha <- 90 - dominant_angle

  cat("Total number of fitted lines (all SED images):", nrow(fit_lines_all), "\n")
  cat("Lines used near dominant direction (±", config$angle_window_deg, "°): ",
      nrow(fit_lines_core), "\n", sep = "")
  cat("Dominant angle of grooves:", round(dominant_angle, 2), "°\n")
  cat("Deviation from 90° (Δα):", round(delta_alpha, 2), "°\n")
  cat("Robust parallelism sigma (MAD-based):", round(sigma_parallel, 2), "°\n")
  # Save parallelism metrics to CSV
  output_dir <- file.path(base_folder, "R_output_dir")
  if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)
  parallel_summary <- data.frame(
    magnification = mag,
    image_path = "all_SED_images",
    n_lines_total = nrow(fit_lines_all),
    n_lines_core = nrow(fit_lines_core),
    dominant_angle_deg = dominant_angle,
    delta_alpha_deg = delta_alpha,
    sigma_parallel_deg = sigma_parallel,
    stringsAsFactors = FALSE
  )
  parallel_file <- file.path(output_dir, sprintf("parallelism_summary_mag%d.csv", mag))
  write.csv(parallel_summary, parallel_file, row.names = FALSE)
  cat("Parallelism metrics saved to:", parallel_file, "\n")
  # 7. Visualization: interactive on screen or saved to PDF
  if (interactive()) {
    par(mfrow = c(1, 2), mar = c(5, 5, 2, 1))
    plot(edges, main = "Edges (custom threshold)")
    plot(img, main = "Fitted lines on edges")
    for (i in seq_along(fit_lines$fit)) {
      fit <- fit_lines$fit[[i]]
      abline(fit, col = "red", lwd = 2)
    }
  } else {
    parallel_plot_file <- file.path(base_folder, "Rplot_parallelism.pdf")
    pdf(parallel_plot_file, width = 8, height = 4)
    par(mfrow = c(1, 2), mar = c(5, 5, 2, 1))
    plot(edges, main = "Edges (custom threshold)")
    plot(img, main = "Fitted lines on edges")
    for (i in seq_along(fit_lines$fit)) {
      fit <- fit_lines$fit[[i]]
      abline(fit, col = "red", lwd = 2)
    }
    dev.off()
    cat("Parallelism plot saved to:", parallel_plot_file, "\n")
  }
} else {
  cat("\n[Parallelism evaluation] No SED images available for parallelism evaluation.\n")
}
}


#----------------------------- distribution_fit_helper ---------------

fit_spacing_distributions <- function(all_mean_distances) {
  all_mean_distances <- as.numeric(all_mean_distances)
  all_mean_distances <- na.omit(all_mean_distances)
  if (length(all_mean_distances) < 5) {
    warning("Not enough data for distribution fitting")
    return(invisible(NULL))
  }

  set.seed(123)
  eps <- sd(all_mean_distances) * 1e-6
  if (is.na(eps) || eps == 0) eps <- 1e-6
  all_mean_distances_jit <- all_mean_distances +
    rnorm(length(all_mean_distances), mean = 0, sd = eps)

  library(MASS)
  library(fitdistrplus)

  fit_normal <- fitdist(all_mean_distances_jit, "norm")
  ks_normal <- suppressWarnings(ks.test(all_mean_distances_jit, "pnorm",
                                       mean = fit_normal$estimate[1],
                                       sd   = fit_normal$estimate[2]))

  fit_exp <- fitdist(all_mean_distances_jit, "exp")
  ks_exp <- suppressWarnings(ks.test(all_mean_distances_jit, "pexp",
                                     rate = 1 / fit_exp$estimate))

  fit_gamma <- fitdist(all_mean_distances_jit, "gamma")
  ks_gamma <- suppressWarnings(ks.test(all_mean_distances_jit, "pgamma",
                                       shape = fit_gamma$estimate[1],
                                       rate  = fit_gamma$estimate[2]))

  fit_lnorm <- fitdist(all_mean_distances_jit, "lnorm")
  ks_lnorm <- suppressWarnings(ks.test(all_mean_distances_jit, "plnorm",
                                       meanlog = fit_lnorm$estimate[1],
                                       sdlog   = fit_lnorm$estimate[2]))

  results <- data.frame(
    Distribution = c("Normal", "Exponential", "Gamma", "Log-Normal"),
    P_value = c(ks_normal$p.value, ks_exp$p.value,
                ks_gamma$p.value, ks_lnorm$p.value)
  )

  print(results)
  best_distribution <- results[which.max(results$P_value), ]
  cat("Best distribution based on p-value:",
      best_distribution$Distribution, "\n")

  hist(all_mean_distances, probability = TRUE, col = "lightblue",
       main = "Empirical Density vs. Fitted Distributions",
      xlab = "Mean Groove Spacing (µm)", border = "white")

  curve(dnorm(x, mean = fit_normal$estimate[1],
              sd = fit_normal$estimate[2]),
        col = "red", lwd = 2, add = TRUE)
  curve(dexp(x, rate = 1 / fit_exp$estimate),
        col = "green", lwd = 2, add = TRUE)
  curve(dgamma(x, shape = fit_gamma$estimate[1],
               rate = fit_gamma$estimate[2]),
        col = "blue", lwd = 2, add = TRUE)
  curve(dlnorm(x, meanlog = fit_lnorm$estimate[1],
               sdlog = fit_lnorm$estimate[2]),
        col = "purple", lwd = 2, add = TRUE)

  legend("topright",
         legend = c("Normal", "Exponential", "Gamma", "Log-Normal"),
         fill   = c("red", "green", "blue", "purple"),
         border = "white")

  invisible(results)
}


run_all_magnifications <- function(mags = c(1000, 750, 500, 400)) {
  for (mag in mags) {
    run_analysis_for_mag(mag)
  }
}


make_summary_plots <- function(mags = c(1000, 750, 500, 400)) {
  output_dir <- file.path(base_folder, "R_output_dir")
  groove_files <- file.path(output_dir,
                            sprintf("groove_summary_mag%d.csv", mags))
  parallel_files <- file.path(output_dir,
                              sprintf("parallelism_summary_mag%d.csv", mags))

  groove_all <- do.call(rbind, lapply(groove_files, read.csv,
                                      stringsAsFactors = FALSE))
  parallel_all <- do.call(rbind, lapply(parallel_files, read.csv,
                                        stringsAsFactors = FALSE))

  summary_all <- merge(groove_all,
                       parallel_all[, c("magnification",
                                        "sigma_parallel_deg")],
                       by = "magnification", all.x = TRUE)

  pdf(file.path(base_folder, "summary_magnification_errors.pdf"),
      width = 7, height = 4)
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

  plot(summary_all$magnification, summary_all$magnification_error_pct,
       type = "b", pch = 19,
      xlab = "Magnification",
      ylab = "Magnification Error [%]",
      main = "Magnification Error vs. Magnification")

  plot(summary_all$magnification, summary_all$groove_spacing_error_um,
       type = "b", pch = 19,
      xlab = "Magnification",
      ylab = "Groove Spacing Error [µm]",
      main = "Groove Spacing Error vs. Magnification")

  dev.off()

  pdf(file.path(base_folder, "summary_parallelism_sigma.pdf"),
      width = 5, height = 4)
  par(mar = c(5, 5, 3, 1))
  plot(summary_all$magnification, summary_all$sigma_parallel_deg,
       type = "b", pch = 19,
      xlab = "Magnification",
      ylab = "Angular Sigma [deg]",
      main = "Groove Parallelism vs. Magnification")
  dev.off()
}
