# From consecutive set of stall sites, identify max position

### WT_fishers ###
orfs <- WT_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- WT_peaks_all[g]$WT_odds
  peaks <- WT_peaks_all[g]$position
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
WT_D4_peaks <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                          peak = final_peaks,
                          WT_D4_odds = final_odds)
WT_D4_peaks[, ID := as.character(base::paste(orf, peak, sep = "_"))]
setkeyv(WT_D4_peaks, c("orf"))
WT_stalling_peaks <- WT_peaks_all[WT_peaks_all$ID %in% WT_D4_peaks$ID]
i <- cbind(match(WT_stalling_peaks$ID, WT_D4_peaks$ID))
WT_stalling_peaks <- cbind(WT_stalling_peaks, peak = WT_D4_peaks[i]$peak)
WT_stalling_peaks <- cbind(WT_stalling_peaks, peak_start = WT_D4_peaks[i]$peak_start)
WT_stalling_peaks <- cbind(WT_stalling_peaks, peak_end = WT_D4_peaks[i]$peak_end)
i <- cbind(match(WT_stalling_peaks$ID, sc_dtA$ID))
WT_stalling_peaks <- cbind(WT_stalling_peaks, motif3 = sc_dtA[i]$motif3)
WT_stalling_peaks <- cbind(WT_stalling_peaks, motif2 = sc_dtA[i]$motif2)
#saveRDS(WT_stalling_peaks, "WT_stalling_peaks.rds")
write.csv(WT_stalling_peaks, "WT_stalling_peaks.csv")

WT_peaks_all_ageDep <- WT_peaks_all[WT_D4_pause > WT_D2_pause]
i <- cbind(match(WT_peaks_all_ageDep$ID, sc_dtA$ID))
WT_peaks_all_ageDep <- cbind(WT_peaks_all_ageDep, motif3 = sc_dtA[i]$motif3)


# WT_stalling_peaks_max <- WT_stalling_peaks[, .SD[which.max(WT_odds)], by = orf]
# write.csv(WT_stalling_peaks_max, "WT_stalling_peaks_max.csv")


### sch9_fishers ###
orfs <- sch9_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- sch9_peaks_all[g]$sch9_odds
  peaks <- sch9_peaks_all[g]$position
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
sch9_D4_peaks <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                            peak = final_peaks,
                            sch9_D4_odds = final_odds)
sch9_D4_peaks[, ID := as.character(base::paste(orf, peak, sep = "_"))]
setkeyv(sch9_D4_peaks, c("orf"))
sch9_stalling_peaks <- sch9_peaks_all[sch9_peaks_all$ID %in% sch9_D4_peaks$ID]
i <- cbind(match(sch9_stalling_peaks$ID, sch9_D4_peaks$ID))
sch9_stalling_peaks <- cbind(sch9_stalling_peaks, peak = sch9_D4_peaks[i]$peak)
sch9_stalling_peaks <- cbind(sch9_stalling_peaks, peak_start = sch9_D4_peaks[i]$peak_start)
sch9_stalling_peaks <- cbind(sch9_stalling_peaks, peak_end = sch9_D4_peaks[i]$peak_end)
i <- cbind(match(sch9_stalling_peaks$ID, sc_dtA$ID))
sch9_stalling_peaks <- cbind(sch9_stalling_peaks, motif3 = sc_dtA[i]$motif3)
sch9_stalling_peaks <- cbind(sch9_stalling_peaks, motif2 = sc_dtA[i]$motif2)
saveRDS(sch9_stalling_peaks, "sch9_stalling_peaks.rds")
write.csv(sch9_stalling_peaks, "sch9_stalling_peaks.csv")
# sch9_stalling_peaks_max <- sch9_stalling_peaks[, .SD[which.max(sch9_odds)], by = orf]
# write.csv(sch9_stalling_peaks_max, "sch9_stalling_peaks_max.csv")

