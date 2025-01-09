library("MALDIquant")
library("MALDIquantForeign")

# Preprocessing (same as your workflow)
data(fiedler2009subset)
spectra <- transformIntensity(fiedler2009subset, method="sqrt")
spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
spectra <- removeBaseline(spectra, method="SNIP", iterations=100)
spectra <- calibrateIntensity(spectra, method="TIC")
spectra <- alignSpectra(spectra, halfWindowSize=20, SNR=2, tolerance=0.002, warpingMethod="lowess")

samples <- factor(sapply(spectra, function(x) metaData(x)$sampleName))
avgSpectra <- averageMassSpectra(spectra, labels=samples, method="mean")
peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=20, SNR=2)
peaks <- binPeaks(peaks, tolerance=0.002)
peaks <- filterPeaks(peaks, minFrequency=0.25)

##################################################################

# Extract aligned spectra data
aligned_mzs <- lapply(avgSpectra, mass)
aligned_intensities <- lapply(avgSpectra, intensity)

# Get detected peak m/z values from the filtered peaks
# peak_mzs <- mass(peaks[[1]])  # Assuming a single avgSpectra; adjust if multiple

peak_mzs <- unique(unlist(lapply(peaks, mass)))
peak_mzs <- sort(peak_mzs)
print(length(peak_mzs))
print(peak_mzs)

# Define bins for each peak
bin_ranges <- lapply(peak_mzs, function(mz) c(mz - 2, mz + 2))
bin_labels <- paste0("bin_", round(peak_mzs, 4))  # Label bins for column names

# Initialize binned data matrix
binned_data <- matrix(0, nrow = length(avgSpectra), ncol = length(bin_ranges))
colnames(binned_data) <- bin_labels

# Sum intensities within bins
for (i in seq_along(avgSpectra)) {
  mzs <- aligned_mzs[[i]]
  intensities <- aligned_intensities[[i]]
  for (j in seq_along(bin_ranges)) {
    range <- bin_ranges[[j]]
    binned_data[i, j] <- sum(intensities[mzs >= range[1] & mzs <= range[2]])
  }
}

# Convert to data frame and export as CSV
binned_data_df <- as.data.frame(binned_data)
row.names(binned_data_df) <- paste0("spectrum_", seq_len(nrow(binned_data_df)))
write.csv(binned_data_df, file = "binned_intensities.csv", row.names = TRUE)

