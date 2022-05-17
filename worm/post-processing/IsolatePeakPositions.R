### Isolate position within stall with max ribosome density ###
orfs <- stalling_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- stalling_peaks_all[g]$odds
  peaks <- stalling_peaks_all[g]$position
  while (length(peaks) > 0) {
    y=peaks[1]
    while (y %in% peaks) {
      peaks_subset1 <- y
      peaks_subset <- c(peaks_subset,peaks_subset1)
      y <- y+1
    } # outputs vector of consecutive positions
    odds_subset <- odds[which(peaks %in% peaks_subset)]
    final_peaks1 <- peaks_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
    final_peaks <- c(final_peaks, final_peaks1) # adds peak to new vector
    final_odds1 <- odds_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
    final_odds <- c(final_odds, final_odds1) # adds peak to new vector
    peak_start1 <- peaks_subset[1] 
    peak_start <- c(peak_start, peak_start1)
    peak_end1 <- peaks_subset[length(peaks_subset)] 
    peak_end <- c(peak_end, peak_end1)
    odds <- odds[which(!peaks %in% peaks_subset)]
    peaks <- peaks[!peaks %in% peaks_subset]
    peaks_subset <- NULL
    odds_subset <- NULL
    gene <- c(gene, g)
  }
}
stalling_peaks1 <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                             peak = final_peaks,
                             odds = final_odds)
stalling_peaks1[, ID := as.character(base::paste(orf, peak, sep = "_"))]
setkeyv(stalling_peaks1, c("orf"))
stalling_peaks <- stalling_peaks_all[stalling_peaks_all$ID %in% stalling_peaks1$ID]
i <- cbind(match(stalling_peaks$ID, stalling_peaks1$ID))
stalling_peaks <- cbind(stalling_peaks, peak = stalling_peaks1[i]$peak)
stalling_peaks <- cbind(stalling_peaks, peak_start = stalling_peaks1[i]$peak_start)
stalling_peaks <- cbind(stalling_peaks, peak_end = stalling_peaks1[i]$peak_end)
#saveRDS(stalling_peaks, "stalling_peaks.rds")
write.csv(stalling_peaks, "stalling_peaks.csv")
# stalling_peaks_max <- stalling_peaks[, .SD[which.max(odds)], by = orf]
# saveRDS(stalling_peaks_max, "stalling_peaks_max.rds")
# write.csv(stalling_peaks_max, "stalling_peaks_max.csv")

