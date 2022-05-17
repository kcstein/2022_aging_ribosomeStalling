### Isolate stall sites using Fisher's exact test ###
his_fishers <- his_dt[WT1_A_rpc >= 1 & WT2_A_rpc >= 1 & 
                        AT1_A_rpc >= 1 & AT2_A_rpc >= 1 & 
                        WT1_A_sum >= 64 & WT2_A_sum >= 64 & 
                        AT1_A_sum >= 64 & AT2_A_sum >= 64 & 
                        position > 20 & stopdist < -20]
his_fishers[, WT := round((WT1_A + WT2_A) / 2)]
his_fishers[, AT := round((AT1_A + AT2_A) / 2)]
his_fishers[, WT_rpc := round((WT1_A_rpc + WT2_A_rpc) / 2)]
his_fishers[, AT_rpc := round((AT1_A_rpc + AT2_A_rpc) / 2)]
his_fishers[, WT_sum := round((WT1_A_sum + WT2_A_sum) / 2)]
his_fishers[, AT_sum := round((AT1_A_sum + AT2_A_sum) / 2)]
his_fishers <- his_fishers[, c(1:7,104:109)]
his_fishers <- his_fishers[!his_fishers$orf %in% his_fishers[(WT_sum - WT) < 0]$orf]
his_fishers <- his_fishers[!his_fishers$orf %in% his_fishers[(AT_sum - AT) < 0]$orf]
View(his_fishers[(WT_sum - WT) < 0]) # Find orfs with position with aberrantly high number of reads
View(his_fishers[(AT_sum - AT) < 0])

for (i in 1:nrow(his_fishers)) {
  print(i)
  counts1 <- matrix(c(his_fishers[i]$AT, his_fishers[i]$WT, 
                      (his_fishers[i]$AT_sum - his_fishers[i]$AT), 
                      (his_fishers[i]$WT_sum - his_fishers[i]$WT)), nrow = 2)
  stalling_test1 <- fisher.test(counts1)
  his_fishers[i, odds := stalling_test1$estimate]
  his_fishers[i, pvalue := stalling_test1$p.value]
}
his_fishers[, padj := p.adjust(pvalue, method = "BH"), by = orf]
#saveRDS(his_fishers, "his_fishers.rds")

peaks_all <- his_fishers[odds > 1 & odds < Inf & padj < 0.05 & 
                          position > 20 & stopdist < -20 &
                          AT > AT_rpc]
#saveRDS(peaks_all, "peaks_all.rds")


orfs <- peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- peaks_all[g]$odds
  peaks <- peaks_all[g]$position
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
AT_peaks <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                           peak = final_peaks,
                           AT_odds = final_odds)
AT_peaks[, ID := as.character(base::paste(orf, peak, sep = "_"))]
setkeyv(AT_peaks, c("orf"))
stalling_peaks <- peaks_all[peaks_all$ID %in% AT_peaks$ID]
i <- cbind(match(stalling_peaks$ID, AT_peaks$ID))
stalling_peaks <- cbind(stalling_peaks, peak = AT_peaks[i]$peak)
stalling_peaks <- cbind(stalling_peaks, peak_start = AT_peaks[i]$peak_start)
stalling_peaks <- cbind(stalling_peaks, peak_end = AT_peaks[i]$peak_end)
#saveRDS(stalling_peaks, "stalling_peaks.rds")


### Codon frequency in stalls ###
aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
proteome_codon <- as.data.table(table(his_dt[position > 20 & stopdist < -20 & residue != "X"]$codon))
i <- cbind(match(stalling_peaks$ID, his_dt$ID))
stalling_peaks <- cbind(stalling_peaks, WT1_A_pause = his_dt[i]$WT1_A_pause)
stalling_peaks <- cbind(stalling_peaks, WT2_A_pause = his_dt[i]$WT2_A_pause)
stalling_peaks <- cbind(stalling_peaks, AT1_A_pause = his_dt[i]$AT1_A_pause)
stalling_peaks <- cbind(stalling_peaks, AT2_A_pause = his_dt[i]$AT2_A_pause)
stalling_peaks[, WT_pause := (WT1_A_pause + WT2_A_pause) / 2]
stalling_peaks[, AT_pause := (AT1_A_pause + AT2_A_pause) / 2]
codon_freq <- as.data.table(table(stalling_peaks[AT_pause > 6]$codon))
codon_freq <- cbind(codon_freq, proteome_codon[, 2])
setnames(codon_freq, c("codon", "stalls", "proteome"))
codon_freq[, stalls_freq := (stalls / sum(stalls))*100]
codon_freq[, proteome_freq := (proteome / sum(proteome))*100]
codon_freq[, div := stalls_freq / proteome_freq]
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(codon_freq$codon, codon_table$codon))
codon_freq <- cbind(codon_freq, residue = codon_table[i]$residue)

codon_freq$codon <- factor(codon_freq$codon, levels = c("GCA", "GCC", "GCG", "GCT",
                                                        "TGC", "TGT",
                                                        "GAC", "GAT",
                                                        "GAA", "GAG",
                                                        "TTC", "TTT",
                                                        "GGA", "GGC", "GGG", "GGT",
                                                        "CAC", "CAT",
                                                        "ATA", "ATC", "ATT",
                                                        "AAA", "AAG",
                                                        "CTA", "CTC", "CTG", "CTT", "TTA", "TTG",
                                                        "ATG",
                                                        "AAC", "AAT",
                                                        "CCA", "CCC", "CCG", "CCT",
                                                        "CAA", "CAG",
                                                        "CGA", "CGC", "CGG", "CGT", "AGA", "AGG",
                                                        "TCA", "TCC", "TCG", "TCT", "AGC", "AGT",
                                                        "ACA", "ACC", "ACG", "ACT",
                                                        "GTA", "GTC", "GTG", "GTT",
                                                        "TGG",
                                                        "TAA", "TAG", "TGA",
                                                        "TAC", "TAT"))

plot <- ggplot(codon_freq[codon != 'TAA' & codon != 'TAG' & codon != 'TGA'], aes(x = codon, y = log2(div), fill = residue)) + geom_col() +
  geom_hline(yintercept = 0, size = 0.2, color = 'black')
plot <- plot + theme_minimal(20, base_line_size = 0.5) + labs(y = "A-site frequency, log2", x = "Codon") +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 5, color = 'black'), axis.ticks.x = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 16, color = "black"), axis.ticks.y = element_line(color = "black"))
ggsave("/Users/KevinStein/Desktop/codonfreq_Green.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Residue frequency in stalls ###
proteome_aa <- table(his_dt[position > 20 & stopdist < -20 & residue != "X"]$residue)
proteome_aa <- as.data.table(proteome_aa[c(1:19,21)])
aa_freq <- table(stalling_peaks[AT_pause > 6]$residue)
aa_freq <- as.data.table(aa_freq[c(1:19,21)])
aa_freq <- cbind(aa_freq, proteome_aa[, 2])
setnames(aa_freq, c("residue", "stalls", "proteome"))
aa_freq[, stalls_freq := (stalls / sum(stalls))*100]
aa_freq[, proteome_freq := (proteome / sum(proteome))*100]
aa_freq[, div := stalls_freq / proteome_freq]

sc_dtA <- readRDS("../../Chronological/doc/sc_dtA.rds")
i <- cbind(match(stalling_peaks$ID, sc_dtA$ID))
stalling_peaks <- cbind(stalling_peaks, motif3 = sc_dtA[i]$motif3)
stalling_peaks[, P := as.factor(substring(motif3,2,2))]
stalling_peaks[, E := as.factor(substring(motif3,1,1))]

aa_freqP <- data.table(table(stalling_peaks[AT_pause > 6]$P))
aa_freqP <- as.data.table(aa_freqP[c(1:19,21)])
aa_freqP <- cbind(aa_freqP, proteome_aa[, 2])
setnames(aa_freqP, c("residue", "stalls", "proteome"))
aa_freqP[, stalls_freq := (stalls / sum(stalls))*100]
aa_freqP[, proteome_freq := (proteome / sum(proteome))*100]
aa_freqP[, div := stalls_freq / proteome_freq]

aa_freqE <- data.table(table(stalling_peaks[AT_pause > 6]$E))
aa_freqE <- as.data.table(aa_freqE[c(1:19,21)])
aa_freqE <- cbind(aa_freqE, proteome_aa[, 2])
setnames(aa_freqE, c("residue", "stalls", "proteome"))
aa_freqE[, stalls_freq := (stalls / sum(stalls))*100]
aa_freqE[, proteome_freq := (proteome / sum(proteome))*100]
aa_freqE[, div := stalls_freq / proteome_freq]


plot <- ggplot(aa_freq[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "A-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/A-siteFreq_Green.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aa_freqP[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "P-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/P-siteFreq_Green.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aa_freqE[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "E-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/E-siteFreq_Green.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

library(Logolas)
bg <- summary(his_fishers[position > 20 & stopdist < -20]$residue) / 
  length(his_fishers[position > 20 & stopdist < -20]$residue)
bg <- bg[c(1:19,21)]

logomaker(stalling_peaks[AT_pause > 6]$motif3, type = "EDLogo", bg = bg, color_seed = 6)

