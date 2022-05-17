ggplot(data = KR6of6_dt[WT_D0_1A_rpc >= 0.1 & WT_D0_2A_rpc >= 0.1 & 
                          WT_D4_1A_rpc >= 0.1 & WT_D4_2A_rpc >= 0.1 & 
                          polybasic > 25 & (polybasic - length) < -25,
                        .(adjusted, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=3, center=T),
                          WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())


aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
residueMean <- sc_dtA[WT_D0_1A_rpc >= 1 & WT_D0_2A_rpc >= 1 & 
                 WT_D2_1A_rpc >= 1 & WT_D2_2A_rpc >= 1 &
                 WT_D4_1A_rpc >= 1 & WT_D4_2A_rpc >= 1 &
                 WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                 WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
                 WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 &
                 WT_D0_1A_coverage > 0.7 & WT_D0_2A_coverage > 0.7 &
                 WT_D4_1A_coverage > 0.7 & WT_D4_2A_coverage > 0.7 & 
                 WT_D2_1A_coverage > 0.7 & WT_D2_2A_coverage > 0.7 &
                 position > 20 & stopdist < -20]

D0_1 <- NULL
D0_2 <- NULL
D4_1 <- NULL
D4_2 <- NULL
p <- NULL
for (i in aa) {
  D0_1a <- mean(residueMean[residue == i]$WT_D0_1A_pause)
  D0_1 <- c(D0_1, D0_1a)
  D0_2a <- mean(residueMean[residue == i]$WT_D0_2A_pause)
  D0_2 <- c(D0_2, D0_2a)
  D4_1a <- mean(residueMean[residue == i]$WT_D4_1A_pause)
  D4_1 <- c(D4_1, D4_1a)
  D4_2a <- mean(residueMean[residue == i]$WT_D4_2A_pause)
  D4_2 <- c(D4_2, D4_2a)
  p2 <- wilcox.test(residueMean[residue == i]$WT_D0_pause, residueMean[residue == i]$WT_D4_pause, alternative = 't')
  p1 <- p2$p.value
  p <- c(p, p1)
}
pauseMean <- data.table(aa, D0_1, D0_2, D4_1, D4_2, p)
pauseMean[, D0 := (D0_1 + D0_2) / 2]
pauseMean[, D4 := (D4_1 + D4_2) / 2]

D0_1 <- NULL
D0_2 <- NULL
D4_1 <- NULL
D4_2 <- NULL
p <- NULL
for (i in codon_table$codon) {
  D0_1a <- mean(residueMean[codon == i]$WT_D0_1A_pause)
  D0_1 <- c(D0_1, D0_1a)
  D0_2a <- mean(residueMean[codon == i]$WT_D0_2A_pause)
  D0_2 <- c(D0_2, D0_2a)
  D4_1a <- mean(residueMean[codon == i]$WT_D4_1A_pause)
  D4_1 <- c(D4_1, D4_1a)
  D4_2a <- mean(residueMean[codon == i]$WT_D4_2A_pause)
  D4_2 <- c(D4_2, D4_2a)
  # p2 <- wilcox.test(residueMean[residue == i]$WT_D0_pause, residueMean[residue == i]$WT_D4_pause, alternative = 't')
  # p1 <- p2$p.value
  # p <- c(p, p1)
}
pauseMeanCodon <- data.table(codon = codon_table$codon, D0_1, D0_2, D4_1, D4_2)
pauseMeanCodon[, D0 := (D0_1 + D0_2) / 2]
pauseMeanCodon[, D4 := (D4_1 + D4_2) / 2]
i <- cbind(match(pauseMeanCodon$codon, codon_table$codon))
pauseMeanCodon <- cbind(pauseMeanCodon, residue = codon_table[i]$residue)


plot <- ggplot(pauseMean, aes(D0, D4, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) + 
  xlim(0.5,1.7) + ylim(0.5,1.7)
plot <- plot + theme_classic(20) + labs(y = "Day 4 pause score", x = "Day 0 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMean, aes(D0, D4)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 3, color = "#E7A427", alpha = 0.6) +
  xlim(0.5,1.7) + ylim(0.5,1.7)
plot <- plot + theme_classic(20) + labs(y = "Day 4 pause score", x = "Day 0 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_yeastPoints.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
cor(pauseMean$D0, pauseMean$D4)

plot <- ggplot(pauseMeanCodon[residue != 'X'], aes(D0, D4, color = residue)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  #geom_point(size = 2.5, alpha = 0.6) +
  geom_text(aes(label = codon), size = 3)
  xlim(0.5,2.5) + ylim(0.5,2.5)
plot <- plot + theme_classic(20) + labs(y = "Day 4 pause score", x = "Day 0 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_yeastPointsCodons.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



plot <- ggplot(pauseMean, aes(D0_1, D0_2, color = aa)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) +
  xlim(0.5,1.7) + ylim(0.5,1.7)
plot <- plot + theme_classic(20) + labs(y = "Day0, rep2 pause score", x = "Day 0, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day0Replicates_yeastResidues.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMean, aes(D4_1, D4_2, color = aa)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) +
  xlim(0.5,1.7) + ylim(0.5,1.7)
plot <- plot + theme_classic(20) + labs(y = "Day4, rep2 pause score", x = "Day 4, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day4Replicates_yeastResidues.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(pauseMeanCodon[residue != 'X'], aes(D0_1, D0_2)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 2.5, alpha = 0.6) +
  xlim(0.5,2.5) + ylim(0.5,2.5)
plot <- plot + theme_classic(20) + labs(y = "Day0, rep2 pause score", x = "Day 0, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day0Replicates_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMeanCodon[residue != 'X'], aes(D4_1, D4_2, color = residue)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  #geom_point(size = 2.5, alpha = 0.6) +
  geom_text(aes(label = codon), size = 3) +
  xlim(0.5,2.5) + ylim(0.5,2.5)
plot <- plot + theme_classic(20) + labs(y = "Day4, rep2 pause score", x = "Day 4, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day4Replicates_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



plot <- ggplot(pause, aes(D0, D4, color = aa)) + geom_text(label = aa, size = 5) + xlim(0.6,1.7) + ylim(0.6,1.7) +
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed")
plot <- plot + theme_classic(20) + labs(y = "Day 4 pause score", x = "Day 0 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

  
plot <- ggplot(pause, aes(D0_1, D0_2, color = aa)) + geom_text(label = aa, size = 5) + xlim(0.6,1.7) + ylim(0.6,1.7) +
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed")
plot <- plot + theme_classic(20) + labs(y = "Rep 1 pause score", x = "Rep 2 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSitesReplicates_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

wilcox.test(temp[residue == "W"]$WT_D0_pause, temp[residue == "W"]$WT_D4_pause, alternative = 'g')
  
