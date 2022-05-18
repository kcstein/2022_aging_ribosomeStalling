aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')

residueMedian <- N2_dtA[D1_1A_rpc >= 1 & D12_1A_rpc >= 1 & 
                 D1_2A_rpc >= 1 & D12_2A_rpc >= 1 & 
                 N2_D6_2A_rpc >= 1 &
                 D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                 D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
                 N2_D6_2A_sum >= 64 & 
                 D1_1A_coverage > 0.7 & D1_2A_coverage > 0.7 &
                 D12_1A_coverage > 0.7 & D12_2A_coverage > 0.7 & 
                 N2_D6_2A_coverage > 0.7 &
                 D1_1A > 0 & D1_2A > 0 &
                 D12_1A > 0 & D12_2A > 0 &
                 N2_D6_2A > 0 & 
                 position > 20 & stopdist < -20]

residueMean <- N2_dtA[D1_1A_rpc >= 1 & D12_1A_rpc >= 1 & 
                          D1_2A_rpc >= 1 & D12_2A_rpc >= 1 & 
                          N2_D6_2A_rpc >= 1 &
                          D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                          D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
                          N2_D6_2A_sum >= 64 & 
                          D1_1A_coverage > 0.7 & D1_2A_coverage > 0.7 &
                          D12_1A_coverage > 0.7 & D12_2A_coverage > 0.7 & 
                          N2_D6_2A_coverage > 0.7 &
                          position > 20 & stopdist < -20]

D1_1 <- NULL
D1_2 <- NULL
D12_1 <- NULL
D12_2 <- NULL
p <- NULL
for (i in aa) {
  D1_1a <- median(residueMedian[residue == i]$D1_1A_pause)
  D1_1 <- c(D1_1, D1_1a)
  D1_2a <- median(residueMedian[residue == i]$D1_2A_pause)
  D1_2 <- c(D1_2, D1_2a)
  D12_1a <- median(residueMedian[residue == i]$D12_1A_pause)
  D12_1 <- c(D12_1, D12_1a)
  D12_2a <- median(residueMedian[residue == i]$D12_2A_pause)
  D12_2 <- c(D12_2, D12_2a)
  p2 <- wilcox.test(residueMedian[residue == i]$D1_pause, residueMedian[residue == i]$D12_pause, alternative = 't')
  p1 <- p2$p.value
  p <- c(p, p1)
}
pauseMedian <- data.table(aa, D1_1, D1_2, D12_1, D12_2, p)
pauseMedian[, D1 := (D1_1 + D1_2) / 2]
pauseMedian[, D12 := (D12_1 + D12_2) / 2]

D1_1 <- NULL
D1_2 <- NULL
D12_1 <- NULL
D12_2 <- NULL
D6 <- NULL
p <- NULL
for (i in aa) {
  D1_1a <- mean(residueMean[residue == i]$D1_1A_pause)
  D1_1 <- c(D1_1, D1_1a)
  D1_2a <- mean(residueMean[residue == i]$D1_2A_pause)
  D1_2 <- c(D1_2, D1_2a)
  D12_1a <- mean(residueMean[residue == i]$D12_1A_pause)
  D12_1 <- c(D12_1, D12_1a)
  D12_2a <- mean(residueMean[residue == i]$D12_2A_pause)
  D12_2 <- c(D12_2, D12_2a)
  D6a <- mean(residueMean[residue == i]$N2_D6_2A_pause)
  D6 <- c(D6, D6a)
  p2 <- wilcox.test(residueMean[residue == i]$D1_pause, residueMean[residue == i]$D12_pause, alternative = 't')
  p1 <- p2$p.value
  p <- c(p, p1)
}
pauseMean <- data.table(aa, D1_1, D1_2, D12_1, D12_2, D6, p)
pauseMean[, D1 := (D1_1 + D1_2) / 2]
pauseMean[, D12 := (D12_1 + D12_2) / 2]

D1_1 <- NULL
D1_2 <- NULL
D12_1 <- NULL
D12_2 <- NULL
p <- NULL
for (i in codon_table$codon) {
  D1_1a <- mean(residueMean[codon == i]$D1_1A_pause)
  D1_1 <- c(D1_1, D1_1a)
  D1_2a <- mean(residueMean[codon == i]$D1_2A_pause)
  D1_2 <- c(D1_2, D1_2a)
  D12_1a <- mean(residueMean[codon == i]$D12_1A_pause)
  D12_1 <- c(D12_1, D12_1a)
  D12_2a <- mean(residueMean[codon == i]$D12_2A_pause)
  D12_2 <- c(D12_2, D12_2a)
  #p2 <- wilcox.test(residueMean[codon == i]$D1_pause, residueMean[codon == i]$D12_pause, alternative = 't')
  #p1 <- p2$p.value
  #p <- c(p, p1)
}
pauseMeanCodon <- data.table(codon = codon_table$codon, D1_1, D1_2, D12_1, D12_2)
pauseMeanCodon[, D1 := (D1_1 + D1_2) / 2]
pauseMeanCodon[, D12 := (D12_1 + D12_2) / 2]
i <- cbind(match(pauseMeanCodon$codon, codon_table$codon))
pauseMeanCodon <- cbind(pauseMeanCodon, residue = codon_table[i]$residue)



plot <- ggplot(pauseMean, aes(D1, D12, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) + 
  xlim(0.85,1.2) + ylim(0.85,1.2)
plot <- plot + theme_classic(20) + labs(y = "Day 12 pause score", x = "Day 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(pauseMeanCodon[residue != 'X'], aes(D1, D12)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 2.5, alpha = 0.8) +
  xlim(0.85,1.25) + ylim(0.85,1.25)
plot <- plot + theme_classic(20) + labs(y = "Day 12 pause score", x = "Day 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_wormsPointsCodons.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(pauseMeanCodon[residue != 'X'], aes(D1_1, D1_2)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 2.5, alpha = 0.6) +
  xlim(0.85,1.25) + ylim(0.85,1.25)
plot <- plot + theme_classic(20) + labs(y = "Day1, rep2 pause score", x = "Day 1, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day1Replicates_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMeanCodon[residue != 'X'], aes(D12_1, D12_2)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 2.5, alpha = 0.6) +
  xlim(0.85,1.25) + ylim(0.85,1.25)
plot <- plot + theme_classic(20) + labs(y = "Day12, rep2 pause score", x = "Day 12, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day12Replicates_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



ggplot(pauseMedian, aes(D1, D12, color = aa)) + geom_text(label = aa, size = 5) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed") +
  #xlim(1.1,1.4) + ylim(1.1,1.4)
  xlim(0.85,1.2) + ylim(0.85,1.2)


plot <- ggplot(pause, aes(D1_1, D1_2, color = aa)) + geom_text(label = aa, size = 5) + xlim(0.9,1.12) + ylim(0.9,1.12) +
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1)
plot <- plot + theme_classic(20) + labs(y = "Rep 1 pause score", x = "Rep 2 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSitesReplicates_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pause, aes(D12_1, D12_2, color = aa)) + geom_text(label = aa, size = 5) + xlim(0.9,1.12) + ylim(0.9,1.12) +
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1)
plot <- plot + theme_classic(20) + labs(y = "Rep 1 pause score", x = "Rep 2 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSitesReplicates_wormsD12.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)