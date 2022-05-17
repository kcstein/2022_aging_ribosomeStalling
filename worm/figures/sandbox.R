### RR
RR <- N2_dtA[motif2 == "RR" & position > 15 & stopdist < -15 &
               D1_1A_rpc >= 1 & D1_2A_rpc >= 1 &
               D12_1A_rpc >= 1 & D12_2A_rpc >= 1 & D1_1A > 0 & D1_2A > 0 & D12_1A > 0 & D12_2A > 0]
RR_coverage <- N2_dtA[motif2 == "RR" & position > 15 & stopdist < -15 &
         D1_1A_rpc >= 1 & D1_2A_rpc >= 1 &
         D12_1A_rpc >= 1 & D12_2A_rpc >= 1 & D1_1A_coverage > 0.5 & D1_2A_coverage > 0.5 & D12_1A_coverage > 0.5 & D12_2A_coverage > 0.5]
ggplot(RR, aes("1", D1_pause)) + geom_boxplot() +
  geom_boxplot(aes("12", D12_pause)) + scale_y_log10()
wilcox.test(RR$D1_pause, RR$D12_pause)
RR1 <- RR[, c(1:2,6)]
setnames(RR1, c("orf", "peak", "name"))
RR_dt <- N2_dtA[RR1, allow.cartesian=TRUE]
RR_dt[, adjusted := position - peak]

### Normalize across polybasic interval ###
test <- RR_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("name"))
setkeyv(RR_dt, c("name"))
test[, D1_1A_rpc_adjusted := mean(D1_1A), by = name]
test[, D1_2A_rpc_adjusted := mean(D1_2A), by = name]
test[, D12_1A_rpc_adjusted := mean(D12_1A), by = name]
test[, D12_2A_rpc_adjusted := mean(D12_2A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(RR_dt$name, test$name))
RR_dt <- cbind(RR_dt, D1_1A_rpc_adjusted = test[i]$D1_1A_rpc_adjusted)
RR_dt <- cbind(RR_dt, D1_2A_rpc_adjusted = test[i]$D1_2A_rpc_adjusted)
RR_dt <- cbind(RR_dt, D12_1A_rpc_adjusted = test[i]$D12_1A_rpc_adjusted)
RR_dt <- cbind(RR_dt, D12_2A_rpc_adjusted = test[i]$D12_2A_rpc_adjusted)
RR_dt[, D1_1A_norm := D1_1A / D1_1A_rpc_adjusted]
RR_dt[, D1_2A_norm := D1_2A / D1_2A_rpc_adjusted]
RR_dt[, D12_1A_norm := D12_1A / D12_1A_rpc_adjusted]
RR_dt[, D12_2A_norm := D12_2A / D12_2A_rpc_adjusted]
# saveRDS(RR_dt, "RR_dt.rds")


plot <- ggplot(data = RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                              D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1,
                            .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=1, center=T),
                              D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7298A') +
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size = 1.25, color = '#E7298A') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size = 1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = "none", legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Worm_RR_smooth1.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                    D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 & 
                    adjusted == 0]$D1_2A_norm,
            RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                    D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 & 
                    adjusted == 0]$D12_2A_norm, alternative = 't')


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


# D1 <- NULL
# D12 <- NULL
# p <- NULL
# for (i in aa) {
#   D1a <- mean(stalling_peaks_ageDep[residue == i]$D1_pause)
#   D1 <- c(D1, D1a)
#   D12a <- mean(stalling_peaks_ageDep[residue == i]$D12_pause)
#   D12 <- c(D12, D12a)
#   p2 <- wilcox.test(stalling_peaks_ageDep[residue == i]$D1_pause, stalling_peaks_ageDep[residue == i]$D12_pause, alternative = 't')
#   p1 <- p2$p.value
#   p <- c(p, p1)
# }
# pause <- data.table(aa, D1, D12, p)
# pause[, ratio := D12/D1]

plot <- ggplot(pauseMean, aes(D1, D12, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) + 
  xlim(0.85,1.2) + ylim(0.85,1.2)
plot <- plot + theme_classic(20) + labs(y = "Day 12 pause score", x = "Day 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMean, aes(D1, D12)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 3, color = "#E7298A", alpha = 0.6) +
  xlim(0.85,1.2) + ylim(0.85,1.2)
plot <- plot + theme_classic(20) + labs(y = "Day 12 pause score", x = "Day 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_wormsPoints.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMeanCodon[residue != 'X'], aes(D1, D12)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 2.5, alpha = 0.8) +
  xlim(0.85,1.25) + ylim(0.85,1.25)
plot <- plot + theme_classic(20) + labs(y = "Day 12 pause score", x = "Day 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_wormsPointsCodons.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(pauseMean, aes(D1_1, D1_2, color = aa)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) +
  xlim(0.85,1.25) + ylim(0.85,1.25)
plot <- plot + theme_classic(20) + labs(y = "Day 1, rep 2 pause score", x = "Day 1, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day1Replicates_wormsResidues.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMean, aes(D12_1, D12_2, color = aa)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) +
  xlim(0.85,1.25) + ylim(0.85,1.25)
plot <- plot + theme_classic(20) + labs(y = "Day12, rep2 pause score", x = "Day 12, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day12Replicates_wormsResidues.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


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

ggplot(N2_fishers, aes(residue, odds)) + geom_boxplot() + coord_cartesian(ylim = c(0,2))

ggplot(stalling_peaks_ageDep_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                                          D6_2A_rpc_adjusted >= 1 & 
                                          D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 &
                                          peak > 25 & (peak - length) < -25,
                                        .(adjusted, D1 = movingAverage(D1_pause, n=2, center=T),
                                          D6 = movingAverage(N2_D6_2A_pause, n=2, center=T),
                                          D12 = movingAverage(D12_pause, n=2, center=T))]) +
  stat_summary(aes(adjusted, D6), fun.y = "mean", geom = "line", size=1.25, color = '#F28CC1') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=1.25, color = '#E7298A') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_color_manual(labels = c("Day 1","Day 6", "Day 12"), values = c("gray40", "#F28CC1", "#E7298A"), name = "") +
  scale_x_continuous(expand = expansion(), limits = c(-25, 25))


temp <- KR6of6_dt[adjusted == 0]
temp1 <- temp[D1_1A > 0 & D1_2A > 0 & D12_1A > 0 & D12_2A > 0]
temp2 <- KR6of6_dt[KR6of6_dt$name %in% temp1$name]

ggplot(data = temp2[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                                  D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1,
                                .(adjusted, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=3, center=T),
                                  D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7298A') +
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size = 1.25, color = '#E7298A') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size = 1.25, color = 'gray40') +
  coord_cartesian(xlim = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = "none", legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Worm_KR6of6.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Correlation of datasets ###

library(corrplot)
library(tidyr)
N2_expression_deseq <- readRDS("doc/N2_expression_deseq.rds")
names(N2_expression_deseq)
M <- cor(N2_expression_deseq[, c(18,20,22,24,26,28,30)])
M <- cor(N2_expression_deseq[, c(18,20,30,22,24)])
head(round(M, 2))
corrplot(M, method="color",  
         type="upper", 
         tl.col="black", tl.srt=45)

pdf(file = "/Users/KevinStein/Desktop/worm_correlation.pdf", width = 6.5, height = 5.5, useDingbats = F)
corrplot(M, method="color",  
         type="upper", 
         tl.col="black", tl.srt=45, cl.align.text = 'l', addCoef.col = "white", number.cex = 0.9
)
dev.off()


### Pausing correlation ###
names(pauseMean)
M <- cor(pauseMean[, c(2,3,6,4,5)])
head(round(M, 2))
corrplot(M, method="color",  
         type="lower", 
         tl.col="black")

pdf(file = "/Users/KevinStein/Desktop/worm_pauseScoreCorrelation.pdf", width = 6.5, height = 5.5, useDingbats = F)
corrplot(M, method="color",  
         type="lower", diag = FALSE,
         tl.col="black", tl.srt=45, cl.align.text = 'l', addCoef.col = "white", number.cex = 1
)
dev.off()


### Example ###
View(N2_dtA[position == 1 & WBgene %in% c("WBGene00010624", "WBGene00017799", "WBGene00020348", "WBGene00004494", "WBGene00004472", "WBGene00004471", "WBGene00001340", "WBGene00004326", "WBGene00004431",
                          "WBGene00016740", "WBGene00004453", "WBGene00007859", "WBGene00004490", "WBGene00004454", "WBGene00004452", "WBGene00004425", "WBGene00004410", "WBGene00004355",
                          "WBGene00017816", "WBGene00006794", "WBGene00017245", "WBGene00010097", "WBGene00004469", "WBGene00000479", "WBGene00004482", "WBGene00000474", "WBGene00004441",
                          "WBGene00015133", "WBGene00012556", "WBGene00004434", "WBGene00015487", "WBGene00016249", "WBGene00004415", "WBGene00004433", "WBGene00004918", "WBGene00004435",
                          "WBGene00004917", "WBGene00015185", "WBGene00004916", "WBGene00012484", "WBGene00016142", "WBGene00004413", "WBGene00015461", "WBGene00016653", "WBGene00016493", "WBGene00022739", "WBGene00004419"), c(70,71)])

plot <- ggplot(data=N2_dtA[WBgene == "WBGene00004472", .(position, D1 = movingAverage(((D1_1A_rpm + D1_2A_rpm) / 2), n=5, center=T), 
                                                 D6 = movingAverage(N2_D6_2A_rpm, n=5, center=T),
                                               D12 = movingAverage(((D12_1A_rpm + D12_2A_rpm) / 2), n=5, center=T))]) + 
  geom_line(aes(position, D1), color = "gray40", size = 1.25) +
  geom_line(aes(position, D6), color = "#F28CC1", size = 1.25) +
  geom_line(aes(position, D12), color = "#E7298A", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/rps-3.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Metagene from start
N2_dtA <- readRDS("doc/N2_dtA.rds")
plot <- ggplot(data = N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & N2_D6_2A_rpc >= 0.5 &
                       D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                       D1_1A_sum >= 64 & D1_2A_sum >= 64 & N2_D6_2A_sum >= 64 &
                       D12_1A_sum >= 64 & D12_2A_sum >= 64 & length >= 300,
                     .(position, D1 = ((D1_1A_pause + D1_2A_pause) / 2),
                       D6 = N2_D6_2A_pause,
                       D12 = ((D12_1A_pause + D12_2A_pause) / 2))]) + xlim(-7, 200) +
  stat_summary(aes(position, D6, color = '#999999'), fun.y = "mean", geom = "line", size=1.25) +
  stat_summary(aes(position, D12, color = '#E31A1C'), fun.y = "mean", geom = "line", size=1.25) +
  stat_summary(aes(position, D1, color = '#1F78B4'), fun.y = "mean", geom = "line", size=1.25) +
  scale_color_manual(labels = c("Day 1","Day 6", "Day 12"), values = c("#1F78B4", "#999999", "#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Norm. ribosome occupancy", x = "Codon position") +
  theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/startMetagene_worm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & N2_D6_2A_rpc >= 0.5 &
                               D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                               D1_1A_sum >= 64 & D1_2A_sum >= 64 & N2_D6_2A_sum >= 64 &
                               D12_1A_sum >= 64 & D12_2A_sum >= 64 & length >= 300,
                             .(position, D1 = ((D1_1A_pause + D1_2A_pause) / 2),
                               D12 = ((D12_1A_pause + D12_2A_pause) / 2))]) + xlim(-7, 200) +
  stat_summary(aes(position, D12, color = '#E7298A'), fun.y = "mean", geom = "line", size=1.25) +
  stat_summary(aes(position, D1, color = 'gray40'), fun.y = "mean", geom = "line", size=1.25) +
  scale_color_manual(labels = c("Day 1", "Day 12"), values = c("gray40", "#E7298A"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Norm. ribosome occupancy", x = "Codon position") +
  theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/startMetagene_worm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Metagene from stop
plot <- ggplot(data = N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & N2_D6_2A_rpc >= 0.5 &
                               D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                               D1_1A_sum >= 64 & D1_2A_sum >= 64 & N2_D6_2A_sum >= 64 &
                               D12_1A_sum >= 64 & D12_2A_sum >= 64 & length >= 300,
                             .(stopdist, D1 = ((D1_1A_pause + D1_2A_pause) / 2),
                               D12 = ((D12_1A_pause + D12_2A_pause) / 2))]) + xlim(-200, 7) +
  stat_summary(aes(stopdist, D1, color = 'gray40'), fun.y = "mean", geom = "line", size=1.25) +
  stat_summary(aes(stopdist, D12, color = '#E7298A'), fun.y = "mean", geom = "line", size=1.25) +
  scale_color_manual(labels = c("Day 1", "Day 12"), values = c("gray40", "#E7298A"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Norm. ribosome occupancy", x = "Codon position") +
  theme(legend.position = c(.3,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stopMetagene_worm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


