# ============================================================
# MALDI-TOF MS preprocessing for a single Bruker fid file
# ============================================================

# -----------------------------
# 1. Load/install packages
# -----------------------------
required_pkgs <- c("MALDIquant", "MALDIquantForeign")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(MALDIquant)
library(MALDIquantForeign)

# -----------------------------
# 2. User settings
# -----------------------------
fid_path <- "/home/stephen/Desktop/Ertapenem_resistance_analysis/ertapenem data/0a0e8488/0_G5/1/1SLin/fid"

output_dir <- "maldi_preprocessing_output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Pre-processing parameters
mass_min <- 2000
mass_max <- 20000
transform_method <- "sqrt"
smooth_method <- "SavitzkyGolay"
half_window_size <- 10
baseline_method <- "SNIP"
baseline_iterations <- 25
calibration_method <- "TIC"
snr_threshold <- 2

# -----------------------------
# 3. Import the fid file
# -----------------------------
# MALDIquantForeign::import can read Bruker flex data directories/files
spectra <- import(fid_path, verbose = FALSE)

if (length(spectra) == 0) {
  stop("No spectra were imported from the provided fid path.")
}

# Use the first spectrum
raw_spectrum <- spectra[[1]]

# -----------------------------
# 4. Trim to 2000-20000 Da
# -----------------------------
trimmed_spectrum <- trim(raw_spectrum, range = c(mass_min, mass_max))

# -----------------------------
# 5. Transformation
# -----------------------------
transformed_spectrum <- transformIntensity(trimmed_spectrum, method = transform_method)

# -----------------------------
# 6. Smoothing
# -----------------------------
smoothed_spectrum <- smoothIntensity(
  transformed_spectrum,
  method = smooth_method,
  halfWindowSize = half_window_size
)

# -----------------------------
# 7. Baseline estimation
# -----------------------------
estimated_baseline <- estimateBaseline(
  smoothed_spectrum,
  method = baseline_method,
  iterations = baseline_iterations
)

# -----------------------------
# 8. Baseline removal
# -----------------------------
baseline_removed_spectrum <- removeBaseline(
  smoothed_spectrum,
  method = baseline_method,
  iterations = baseline_iterations
)

# -----------------------------
# 9. Intensity calibration (TIC)
# -----------------------------
calibrated_spectrum <- calibrateIntensity(
  baseline_removed_spectrum,
  method = calibration_method
)

# -----------------------------
# 10. Peak detection
# -----------------------------
peaks <- detectPeaks(
  calibrated_spectrum,
  SNR = snr_threshold,
  halfWindowSize = half_window_size
)

# -----------------------------
# 11. Save detected peaks table
# -----------------------------
peak_df <- data.frame(
  mz = mass(peaks),
  intensity = intensity(peaks)
)

write.csv(
  peak_df,
  file = file.path(output_dir, "detected_peaks.csv"),
  row.names = FALSE
)

# -----------------------------
# 12. Plotting helper
# -----------------------------
save_png <- function(filename, width = 1800, height = 1200, res = 200, expr) {
  png(filename = filename, width = width, height = height, res = res)
  expr
  dev.off()
}

# Also save a single multi-page PDF
pdf(file.path(output_dir, "preprocessing_plots.pdf"), width = 10, height = 6)

# -----------------------------
# 13. Plot 1: Raw trimmed spectrum
# -----------------------------
plot(
  trimmed_spectrum,
  main = "Raw Spectrum (Trimmed to 2000-20000 Da)",
  xlab = "m/z (Da)",
  ylab = "Intensity"
)

save_png(file.path(output_dir, "01_raw_trimmed_spectrum.png"), expr = {
  plot(
    trimmed_spectrum,
    main = "Raw Spectrum (Trimmed to 2000-20000 Da)",
    xlab = "m/z (Da)",
    ylab = "Intensity"
  )
})

# -----------------------------
# 14. Plot 2: After transformation and smoothing
# -----------------------------
plot(
  transformed_spectrum,
  main = "After Transformation and Smoothing",
  xlab = "m/z (Da)",
  ylab = "Intensity",
  col = "grey60"
)
lines(smoothed_spectrum, col = "blue", lwd = 1.5)
legend(
  "topright",
  legend = c("Transformed", "Smoothed"),
  col = c("grey60", "blue"),
  lwd = c(1, 1.5),
  bty = "n"
)

save_png(file.path(output_dir, "02_transformed_and_smoothed.png"), expr = {
  plot(
    transformed_spectrum,
    main = "After Transformation and Smoothing",
    xlab = "m/z (Da)",
    ylab = "Intensity",
    col = "grey60"
  )
  lines(smoothed_spectrum, col = "blue", lwd = 1.5)
  legend(
    "topright",
    legend = c("Transformed", "Smoothed"),
    col = c("grey60", "blue"),
    lwd = c(1, 1.5),
    bty = "n"
  )
})

# -----------------------------
# 15. Plot 3: Baseline estimate over spectrum
# -----------------------------
plot(
  smoothed_spectrum,
  main = paste0("Baseline Estimation (", baseline_method, ", iterations=", baseline_iterations, ")"),
  xlab = "m/z (Da)",
  ylab = "Intensity"
)
lines(estimated_baseline, col = "red", lwd = 2)
legend(
  "topright",
  legend = c("Smoothed spectrum", "Estimated baseline"),
  col = c("black", "red"),
  lwd = c(1, 2),
  bty = "n"
)

save_png(file.path(output_dir, "03_baseline_estimated.png"), expr = {
  plot(
    smoothed_spectrum,
    main = paste0("Baseline Estimation (", baseline_method, ", iterations=", baseline_iterations, ")"),
    xlab = "m/z (Da)",
    ylab = "Intensity"
  )
  lines(estimated_baseline, col = "red", lwd = 2)
  legend(
    "topright",
    legend = c("Smoothed spectrum", "Estimated baseline"),
    col = c("black", "red"),
    lwd = c(1, 2),
    bty = "n"
  )
})

# -----------------------------
# 16. Plot 4: Baseline removed
# -----------------------------
plot(
  baseline_removed_spectrum,
  main = "After Baseline Removal",
  xlab = "m/z (Da)",
  ylab = "Intensity"
)

save_png(file.path(output_dir, "04_baseline_removed.png"), expr = {
  plot(
    baseline_removed_spectrum,
    main = "After Baseline Removal",
    xlab = "m/z (Da)",
    ylab = "Intensity"
  )
})

# -----------------------------
# 17. Plot 5: After TIC calibration
# -----------------------------
plot(
  calibrated_spectrum,
  main = "After TIC Intensity Calibration",
  xlab = "m/z (Da)",
  ylab = "Calibrated intensity"
)

save_png(file.path(output_dir, "05_tic_calibrated.png"), expr = {
  plot(
    calibrated_spectrum,
    main = "After TIC Intensity Calibration",
    xlab = "m/z (Da)",
    ylab = "Calibrated intensity"
  )
})

# -----------------------------
# 18. Plot 6: Peak detection
# -----------------------------
plot(
  calibrated_spectrum,
  main = paste0("Peak Detection (SNR=", snr_threshold, ", halfWindowSize=", half_window_size, ")"),
  xlab = "m/z (Da)",
  ylab = "Calibrated intensity"
)
points(peaks, col = "red", pch = 4, cex = 0.8)

save_png(file.path(output_dir, "06_peak_detection.png"), expr = {
  plot(
    calibrated_spectrum,
    main = paste0("Peak Detection (SNR=", snr_threshold, ", halfWindowSize=", half_window_size, ")"),
    xlab = "m/z (Da)",
    ylab = "Calibrated intensity"
  )
  points(peaks, col = "red", pch = 4, cex = 0.8)
})

dev.off()

# -----------------------------
# 19. Print summary
# -----------------------------
cat("Preprocessing complete.\n")
cat("Output directory:", normalizePath(output_dir), "\n")
cat("Number of detected peaks:", length(peaks), "\n")
cat("Detected peaks saved to:", file.path(output_dir, "detected_peaks.csv"), "\n")
cat("Plots saved as PNGs and as preprocessing_plots.pdf\n")