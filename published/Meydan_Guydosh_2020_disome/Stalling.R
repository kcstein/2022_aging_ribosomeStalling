library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
disome_dt <- readRDS("disome_dt.rds")

### Identify max disome position in each gene ###
disome_positions <- disome_dt[, .SD[which.max(WT_di)], by = orf]
disome_positions[, ratio := WT_di / WT_mono]
disome_positions <- disome_positions[ratio > 2 & WT_di < Inf & ratio < Inf]
write.csv(disome_positions[, c(1,2,6,46)], "/Users/KevinStein/Desktop/disome_positions.csv")

### Codon frequency in stalls ###
aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
#mm_dt <- readRDS("mm_dt.rds")
proteome_codon <- as.data.table(table(mm_dt[position > 20 & stopdist < -20 & residue != "X"]$codon))
#saveRDS(proteome_codon, "proteome_codon.rds")
Kstalling_peaks <- readRDS("Kstalling_peaks.rds")
Lstalling_peaks <- readRDS("Lstalling_peaks.rds")
proteome_codon <- readRDS("proteome_codon.rds")
Kcodon_freq <- as.data.table(table(Kstalling_peaks[K32m_pause > 4]$codon))
Kcodon_freq <- cbind(Kcodon_freq, proteome_codon[, 2])
setnames(Kcodon_freq, c("codon", "stalls", "proteome"))
Kcodon_freq[, stalls_freq := (stalls / sum(stalls))*100]
Kcodon_freq[, proteome_freq := (proteome / sum(proteome))*100]
Kcodon_freq[, div := stalls_freq / proteome_freq]
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(Kcodon_freq$codon, codon_table$codon))
Kcodon_freq <- cbind(Kcodon_freq, residue = codon_table[i]$residue)

Kcodon_freq$codon <- factor(Kcodon_freq$codon, levels = c("GCA", "GCC", "GCG", "GCT",
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

Lcodon_freq <- as.data.table(table(Lstalling_peaks[L32m_pause > 4]$codon))
Lcodon_freq <- cbind(Lcodon_freq, proteome_codon[, 2])
setnames(Lcodon_freq, c("codon", "stalls", "proteome"))
Lcodon_freq[, stalls_freq := (stalls / sum(stalls))*100]
Lcodon_freq[, proteome_freq := (proteome / sum(proteome))*100]
Lcodon_freq[, div := stalls_freq / proteome_freq]
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(Lcodon_freq$codon, codon_table$codon))
Lcodon_freq <- cbind(Lcodon_freq, residue = codon_table[i]$residue)

Lcodon_freq$codon <- factor(Lcodon_freq$codon, levels = c("GCA", "GCC", "GCG", "GCT",
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

plot <- ggplot(Kcodon_freq[codon != 'TAA' & codon != 'TAG' & codon != 'TGA'], aes(x = codon, y = log2(div), fill = residue)) + geom_col() +
  geom_hline(yintercept = 0, size = 0.2, color = 'black')
plot <- plot + theme_minimal(20, base_line_size = 0.5) + labs(y = "A-site frequency, log2", x = "Codon") +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 5, color = 'black'), axis.ticks.x = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 16, color = "black"), axis.ticks.y = element_line(color = "black"))
ggsave("/Users/KevinStein/Desktop/codonfreq_mouseKidney.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(Lcodon_freq[codon != 'TAA' & codon != 'TAG' & codon != 'TGA'], aes(x = codon, y = log2(div), fill = residue)) + geom_col() +
  geom_hline(yintercept = 0, size = 0.2, color = 'black')
plot <- plot + theme_minimal(20, base_line_size = 0.5) + labs(y = "A-site frequency, log2", x = "Codon") +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 5, color = 'black'), axis.ticks.x = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 16, color = "black"), axis.ticks.y = element_line(color = "black"))
ggsave("/Users/KevinStein/Desktop/codonfreq_mouseLiver.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Residue frequency in stalls ###
#rn_dt <- readRDS("doc/rn_dt.rds")
#proteome_aa <- table(rn_dt[position > 20 & stopdist < -20 & residue != "X"]$residue)
#proteome_aa <- as.data.table(proteome_aa[c(1:19,21)])
#saveRDS(proteome_aa, "proteome_aa.rds")
proteome_aa <- readRDS("proteome_aa.rds")
Kaa_freq <- table(Kstalling_peaks$residue)
Kaa_freq <- as.data.table(Kaa_freq[c(1:19,21)])
Kaa_freq <- cbind(Kaa_freq, proteome_aa[, 2])
setnames(Kaa_freq, c("residue", "stalls", "proteome"))
Kaa_freq[, stalls_freq := (stalls / sum(stalls))*100]
Kaa_freq[, proteome_freq := (proteome / sum(proteome))*100]
Kaa_freq[, div := stalls_freq / proteome_freq]

Kstalling_peaks[, P := as.factor(substring(motif3,2,2))]
Kstalling_peaks[, E := as.factor(substring(motif3,1,1))]

Kaa_freqP <- data.table(table(Kstalling_peaks[K32m_pause > 4]$P))
Kaa_freqP <- cbind(Kaa_freqP, proteome_aa[, 2])
setnames(Kaa_freqP, c("residue", "stalls", "proteome"))
Kaa_freqP[, stalls_freq := (stalls / sum(stalls))*100]
Kaa_freqP[, proteome_freq := (proteome / sum(proteome))*100]
Kaa_freqP[, div := stalls_freq / proteome_freq]

Kaa_freqE <- data.table(table(Kstalling_peaks[K32m_pause > 4]$E))
Kaa_freqE <- as.data.table(Kaa_freqE[c(1:19,21)])
Kaa_freqE <- cbind(Kaa_freqE, proteome_aa[, 2])
setnames(Kaa_freqE, c("residue", "stalls", "proteome"))
Kaa_freqE[, stalls_freq := (stalls / sum(stalls))*100]
Kaa_freqE[, proteome_freq := (proteome / sum(proteome))*100]
Kaa_freqE[, div := stalls_freq / proteome_freq]


Laa_freq <- table(Lstalling_peaks[L32m_pause > 4]$residue)
Laa_freq <- as.data.table(Laa_freq[c(1:19,21)])
Laa_freq <- cbind(Laa_freq, proteome_aa[, 2])
setnames(Laa_freq, c("residue", "stalls", "proteome"))
Laa_freq[, stalls_freq := (stalls / sum(stalls))*100]
Laa_freq[, proteome_freq := (proteome / sum(proteome))*100]
Laa_freq[, div := stalls_freq / proteome_freq]

Lstalling_peaks[, P := as.factor(substring(motif3,2,2))]
Lstalling_peaks[, E := as.factor(substring(motif3,1,1))]

Laa_freqP <- data.table(table(Lstalling_peaks[L32m_pause > 4]$P))
Laa_freqP <- cbind(Laa_freqP, proteome_aa[, 2])
setnames(Laa_freqP, c("residue", "stalls", "proteome"))
Laa_freqP[, stalls_freq := (stalls / sum(stalls))*100]
Laa_freqP[, proteome_freq := (proteome / sum(proteome))*100]
Laa_freqP[, div := stalls_freq / proteome_freq]

Laa_freqE <- data.table(table(Lstalling_peaks[L32m_pause > 4]$E))
Laa_freqE <- cbind(Laa_freqE, proteome_aa[, 2])
setnames(Laa_freqE, c("residue", "stalls", "proteome"))
Laa_freqE[, stalls_freq := (stalls / sum(stalls))*100]
Laa_freqE[, proteome_freq := (proteome / sum(proteome))*100]
Laa_freqE[, div := stalls_freq / proteome_freq]


plot <- ggplot(Kaa_freq[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "A-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/A-siteFreq_mouseKidney.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(Kaa_freqP[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "P-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/P-siteFreq_mouseKidney.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(Kaa_freqE[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "E-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/E-siteFreq_mouseKidney.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(Laa_freq[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "A-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/A-siteFreq_mouseLiver.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(Laa_freqP[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "P-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/P-siteFreq_mouseLiver.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(Laa_freqE[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "E-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/E-siteFreq_mouseLiver.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Average pausing by residue
rn_dt <- readRDS("doc/rn_dt.rds")
residueMean <- rn_dt[K3m_1_A_rpc >= 1 & K3m_2_A_rpc >= 1 & K3m_3_A_rpc >= 1 & 
                       K32m_1_A_rpc >= 1 & K32m_2_A_rpc >= 1 & K32m_3_A_rpc >= 1 &
                       K3m_1_A_sum >= 64 & K3m_2_A_sum >= 64 & K3m_3_A_sum >= 64 &
                       K32m_1_A_sum >= 64 & K32m_2_A_sum >= 64 & K32m_3_A_sum >= 64 &
                       position > 20 & stopdist < -20]
residueMedian <- rn_dt[K3m_1_A > 0 & K3m_2_A > 0 & K3m_3_A > 0 & 
                         K32m_1_A > 0 & K32m_2_A > 0 & K32m_3_A > 0 &
                         position > 20 & stopdist < -20]
residueMean[, K3m_pause := (K3m_1_A_pause + K3m_2_A_pause + K3m_3_A_pause) / 3]
residueMean[, K32m_pause := (K32m_1_A_pause + K32m_2_A_pause + K32m_3_A_pause) / 3]
residueMedian[, K3m_pause := (K3m_1_A_pause + K3m_2_A_pause + K3m_3_A_pause) / 3]
residueMedian[, K32m_pause := (K32m_1_A_pause + K32m_2_A_pause + K32m_3_A_pause) / 3]


K3m_1 <- NULL
K3m_2 <- NULL
K3m_3 <- NULL
K32m_1 <- NULL
K32m_2 <- NULL
K32m_3 <- NULL
p <- NULL
for (i in aa) {
  K3m_1a <- mean(residueMean[residue == i]$K3m_1_A_pause)
  K3m_1 <- c(K3m_1, K3m_1a)
  K3m_2a <- mean(residueMean[residue == i]$K3m_2_A_pause)
  K3m_2 <- c(K3m_2, K3m_2a)
  K3m_3a <- mean(residueMean[residue == i]$K3m_3_A_pause)
  K3m_3 <- c(K3m_3, K3m_3a)
  K32m_1a <- mean(residueMean[residue == i]$K32m_1_A_pause)
  K32m_1 <- c(K32m_1, K32m_1a)
  K32m_2a <- mean(residueMean[residue == i]$K32m_2_A_pause)
  K32m_2 <- c(K32m_2, K32m_2a)
  K32m_3a <- mean(residueMean[residue == i]$K32m_3_A_pause)
  K32m_3 <- c(K32m_3, K32m_3a)
  p2 <- wilcox.test(residueMean[residue == i]$K3m_pause, residueMean[residue == i]$K32m_pause, alternative = 't')
  p1 <- p2$p.value
  p <- c(p, p1)
}
pauseMean <- data.table(aa, K3m_1, K3m_2, K3m_3, K32m_1, K32m_2, K32m_3, p)
pauseMean[, K3m := (K3m_1 + K3m_2 + K3m_3) / 3]
pauseMean[, K32m := (K32m_1 + K32m_2 + K32m_3) / 3]

plot <- ggplot(pauseMean, aes(K3m, K32m, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) + 
  xlim(0.5,2.2) + ylim(0.5,2.2)
plot <- plot + theme_classic(20) + labs(y = "K32m pause score", x = "K3m pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_mouseKidney.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


K3m_1 <- NULL
K3m_2 <- NULL
K3m_3 <- NULL
K32m_1 <- NULL
K32m_2 <- NULL
K32m_3 <- NULL
p <- NULL
for (i in aa) {
  K3m_1a <- median(na.omit(residueMedian[residue == i]$K3m_1_A_pause))
  K3m_1 <- c(K3m_1, K3m_1a)
  K3m_2a <- median(na.omit(residueMedian[residue == i]$K3m_2_A_pause))
  K3m_2 <- c(K3m_2, K3m_2a)
  K3m_3a <- median(na.omit(residueMedian[residue == i]$K3m_3_A_pause))
  K3m_3 <- c(K3m_3, K3m_3a)
  K32m_1a <- median(na.omit(residueMedian[residue == i]$K32m_1_A_pause))
  K32m_1 <- c(K32m_1, K32m_1a)
  K32m_2a <- median(na.omit(residueMedian[residue == i]$K32m_2_A_pause))
  K32m_2 <- c(K32m_2, K32m_2a)
  K32m_3a <- median(na.omit(residueMedian[residue == i]$K32m_3_A_pause))
  K32m_3 <- c(K32m_3, K32m_3a)
  p2 <- wilcox.test(residueMedian[residue == i]$K3m_pause, residueMedian[residue == i]$K32m_pause, alternative = 't')
  p1 <- p2$p.value
  p <- c(p, p1)
}
pauseMedian <- data.table(aa, K3m_1, K3m_2, K3m_3, K32m_1, K32m_2, K32m_3, p)
pauseMedian[, K3m := (K3m_1 + K3m_2 + K3m_3) / 3]
pauseMedian[, K32m := (K32m_1 + K32m_2 + K32m_3) / 3]

plot <- ggplot(pauseMedian, aes(K3m, K32m, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) + 
  xlim(0.5,2.1) + ylim(0.5,2.1)
plot <- plot + theme_classic(20) + labs(y = "K32m pause score", x = "K3m pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_mouseKidney_median.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
