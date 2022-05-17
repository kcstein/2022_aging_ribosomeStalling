library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')


### Examine genome-wide pausing ###
sc_dtA <- readRDS("doc/sc_dtA.rds")
ggplot(sc_dtA[WT_D0_1A_rpc >= 1 & WT_D4_1A_rpc >= 1 & residue != "X"], aes(residue, WT_D0_1A_pause)) + geom_violin(aes(residue, WT_D4_1A_pause), color = "red", fill = NA) + geom_violin(fill = NA) +
  scale_y_log10() + coord_cartesian(ylim = c(0.01,400))
ggplot(sc_dtA[WT_D0_2A_rpc >= 1 & WT_D4_2A_rpc >= 1 & residue != "X"], aes(residue, WT_D0_2A_pause)) + geom_violin(aes(residue, WT_D4_2A_pause), color = "red", fill = NA) + geom_violin(fill = NA) +
  scale_y_log10() + coord_cartesian(ylim = c(0.01,400))

ggplot(sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 & 
                WT_D2_1A_rpc >= 0.5 & WT_D2_2A_rpc >= 0.5 &
                WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 &
                WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
                WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 &
                position > 100 & stopdist < -30], 
       aes(WT_D0_1A_pause)) + stat_ecdf(geom = "step", color = "black") +
  scale_x_log10() + coord_cartesian(xlim = c(0.01,20)) +
  stat_ecdf(aes(WT_D0_2A_pause), geom = "step", color = "black") +
  stat_ecdf(aes(WT_D2_1A_pause), geom = "step", color = "blue") +
  stat_ecdf(aes(WT_D2_2A_pause), geom = "step", color = "orange") + # Rep 2 is better
  stat_ecdf(aes(WT_D4_1A_pause), geom = "step", color = "red") +
  stat_ecdf(aes(WT_D4_2A_pause), geom = "step", color = "red")

plot <- ggplot(sc_dtA[WT_D0_1A_rpc >= 1 & WT_D0_2A_rpc >= 1 & 
                        WT_D2_1A_rpc >= 1 & WT_D2_2A_rpc >= 1 &
                        WT_D4_1A_rpc >= 1 & WT_D4_2A_rpc >= 1 &
                        WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                        WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
                        WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 &
                 position > 100 & stopdist < -30], 
       aes(WT_D0_pause)) + stat_ecdf(geom = "step", color = "black") +
  scale_x_log10() + coord_cartesian(xlim = c(0.01,20)) +
  stat_ecdf(aes(WT_D2_pause), geom = "step", color = "blue") +
  stat_ecdf(aes(WT_D4_pause), geom = "step", color = "red")
plot <- plot + theme_classic(16) + labs(y = "Cumulative fraction", x = "Pause score") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/globalPausing_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


ggplot(sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 &
                WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 &
                WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & position > 30 & stopdist < -5], 
       aes(WT_D0_pause)) + stat_ecdf(geom = "step", color = "black", size = 1.25) +
  stat_ecdf(aes(WT_D4_pause), geom = "step", color = "red", size = 1.25) +
  scale_x_log10() + coord_cartesian(xlim = c(0.05,30))

ggplot(sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 &
                WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 &
                WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & position > 30 & stopdist < -5 & residue == "R"], 
       aes(WT_D0_pause)) + stat_ecdf(geom = "step", color = "black", size = 1.25) +
  stat_ecdf(aes(WT_D4_pause), geom = "step", color = "red", size = 1.25) +
  scale_x_log10() + coord_cartesian(xlim = c(0.05,30))


### Identify stall sites using Fisher's exact test ###
sc_fishers <- sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 &
                       WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 &
                       WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                       WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64]
sc_fishers[, WT_D0 := round((WT_D0_1A + WT_D0_2A) / 2)]
sc_fishers[, WT_D4 := round((WT_D4_1A + WT_D4_2A) / 2)]
sc_fishers[, WT_D0_rpc := round((WT_D0_1A_rpc + WT_D0_2A_rpc) / 2)]
sc_fishers[, WT_D4_rpc := round((WT_D4_1A_rpc + WT_D4_2A_rpc) / 2)]
sc_fishers[, WT_D0_sum := round((WT_D0_1A_sum + WT_D0_2A_sum) / 2)]
sc_fishers[, WT_D4_sum := round((WT_D4_1A_sum + WT_D4_2A_sum) / 2)]
sc_fishers[, sch9_D0 := round((sch9_D0_1A + sch9_D0_2A) / 2)]
sc_fishers[, sch9_D4 := round((sch9_D4_1A + sch9_D4_2A) / 2)]
sc_fishers[, sch9_D0_rpc := round((sch9_D0_1A_rpc + sch9_D0_2A_rpc) / 2)]
sc_fishers[, sch9_D4_rpc := round((sch9_D4_1A_rpc + sch9_D4_2A_rpc) / 2)]
sc_fishers[, sch9_D0_sum := round((sch9_D0_1A_sum + sch9_D0_2A_sum) / 2)]
sc_fishers[, sch9_D4_sum := round((sch9_D4_1A_sum + sch9_D4_2A_sum) / 2)]
sc_fishers <- sc_fishers[, c(1:7,93:96,98:109)]
sc_fishers <- sc_fishers[!sc_fishers$orf %in% sc_fishers[(WT_D0_sum - WT_D0) < 0]$orf]
sc_fishers <- sc_fishers[!sc_fishers$orf %in% sc_fishers[(WT_D4_sum - WT_D4) < 0]$orf]
sc_fishers <- sc_fishers[!sc_fishers$orf %in% sc_fishers[(sch9_D0_sum - sch9_D0) < 0]$orf]
sc_fishers <- sc_fishers[!sc_fishers$orf %in% sc_fishers[(sch9_D4_sum - sch9_D4) < 0]$orf]
View(sc_fishers[(WT_D0_sum - WT_D0) < 0]) # Find orfs with position with aberrantly high number of reads
View(sc_fishers[(WT_D4_sum - WT_D4) < 0])
View(sc_fishers[(sch9_D0_sum - sch9_D0) < 0])
View(sc_fishers[(sch9_D4_sum - sch9_D4) < 0])

for (i in 1:nrow(sc_fishers)) {
  print(i)
  counts1 <- matrix(c(sc_fishers[i]$WT_D4, sc_fishers[i]$WT_D0, 
                      (sc_fishers[i]$WT_D4_sum - sc_fishers[i]$WT_D4), 
                      (sc_fishers[i]$WT_D0_sum - sc_fishers[i]$WT_D0)), nrow = 2)
  counts2 <- matrix(c(sc_fishers[i]$sch9_D4, sc_fishers[i]$sch9_D0, 
                      (sc_fishers[i]$sch9_D4_sum - sc_fishers[i]$sch9_D4), 
                      (sc_fishers[i]$sch9_D0_sum - sc_fishers[i]$sch9_D0)), nrow = 2)
  stalling_test1 <- fisher.test(counts1)
  stalling_test2 <- fisher.test(counts2)
  sc_fishers[i, WT_odds := stalling_test1$estimate]
  sc_fishers[i, WT_pvalue := stalling_test1$p.value]
  sc_fishers[i, sch9_odds := stalling_test2$estimate]
  sc_fishers[i, sch9_pvalue := stalling_test2$p.value]
}
sc_fishers[, WT_padj := p.adjust(WT_pvalue, method = "BH"), by = orf]
sc_fishers[, sch9_padj := p.adjust(sch9_pvalue, method = "BH"), by = orf]
i <- cbind(match(sc_fishers$ID, sc_dtA$ID))
sc_fishers <- cbind(sc_fishers, WT_D2_1A = sc_dtA[i]$WT_D2_1A)
sc_fishers <- cbind(sc_fishers, WT_D2_2A = sc_dtA[i]$WT_D2_2A)
sc_fishers <- cbind(sc_fishers, WT_D2_1A_sum = sc_dtA[i]$WT_D2_1A_sum)
sc_fishers <- cbind(sc_fishers, WT_D2_2A_sum = sc_dtA[i]$WT_D2_2A_sum)
sc_fishers <- cbind(sc_fishers, WT_D2_1A_rpc = sc_dtA[i]$WT_D2_1A_rpc)
sc_fishers <- cbind(sc_fishers, WT_D2_2A_rpc = sc_dtA[i]$WT_D2_2A_rpc)
sc_fishers <- cbind(sc_fishers, WT_D2_1A_pause = sc_dtA[i]$WT_D2_1A_pause)
sc_fishers <- cbind(sc_fishers, WT_D2_2A_pause = sc_dtA[i]$WT_D2_2A_pause)
sc_fishers <- cbind(sc_fishers, sch9_D2_1A = sc_dtA[i]$sch9_D2_1A)
sc_fishers <- cbind(sc_fishers, sch9_D2_2A = sc_dtA[i]$sch9_D2_2A)
sc_fishers <- cbind(sc_fishers, sch9_D2_1A_sum = sc_dtA[i]$sch9_D2_1A_sum)
sc_fishers <- cbind(sc_fishers, sch9_D2_2A_sum = sc_dtA[i]$sch9_D2_2A_sum)
sc_fishers <- cbind(sc_fishers, sch9_D2_1A_rpc = sc_dtA[i]$sch9_D2_1A_rpc)
sc_fishers <- cbind(sc_fishers, sch9_D2_2A_rpc = sc_dtA[i]$sch9_D2_2A_rpc)
sc_fishers <- cbind(sc_fishers, sch9_D2_1A_pause = sc_dtA[i]$sch9_D2_1A_pause)
sc_fishers <- cbind(sc_fishers, sch9_D2_2A_pause = sc_dtA[i]$sch9_D2_2A_pause)
sc_fishers <- cbind(sc_fishers, WT_D2_pause = sc_dtA[i]$WT_D2_pause)
sc_fishers <- cbind(sc_fishers, sch9_D2_pause = sc_dtA[i]$sch9_D2_pause)
#saveRDS(sc_fishers, "sc_fishers.rds")

ggplot(WT_stalling_peaks, aes(WT_D4_pause, sch9_D4_pause)) + geom_point() + 
  scale_x_log10(limits = c(0.1,20)) + scale_y_log10(limits = c(0.1,20)) +
  geom_abline(slope = 1, intercept = 0, color = "red")
ggplot(WT_stalling_peaks, aes("WT", WT_D4_pause)) + geom_boxplot() +
  geom_boxplot(aes("sch9", sch9_D4_pause)) + coord_cartesian(ylim = c(0,5))
scale_y_log10(limits = c(0.1,20))
test <- wilcox.test(WT_stalling_peaks$WT_D4_pause, WT_stalling_peaks$sch9_D4_pause, alternative = 't')
test$p.value

WT_stalling_peaks <- readRDS("WT_stalling_peaks.rds")
plot <- ggplot(WT_stalling_peaks[sch9_odds > 0], aes("1WT", WT_odds)) + geom_boxplot(notch = T, fill = "#E7A427") +
  geom_boxplot(aes("2sch9", sch9_odds), fill = "#00CED1", notch = T) + 
  coord_capped_cart(left = "both", ylim = c(0,5))
plot <- plot + theme_classic(20) + labs(y = "odds ratio (Day 4 / Day 0)", x = "") +
  theme(axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Day4_WTvsSch9dark.pdf", plot, width = 3, height = 6, dpi = 300, useDingbats = F)

test <- wilcox.test(WT_stalling_peaks[sch9_D4_pause > 0]$WT_odds, WT_stalling_peaks$sch9_odds, alternative = 't')
test$p.value


sc_fishers <- readRDS("sc_fishers.rds")
WT_peaks_all <- sc_fishers[WT_odds > 1 & WT_odds < Inf & WT_padj < 0.05 & 
                             position > 20 & stopdist < -20 &
                             WT_D4 > WT_D4_rpc]
#saveRDS(WT_peaks_all, "WT_peaks_all.rds")
sch9_peaks_all <- sc_fishers[sch9_odds > 1 & sch9_odds < Inf & sch9_padj < 0.05 & 
                               position > 20 & stopdist < -20 &
                               sch9_D4 > sch9_D4_rpc]
#saveRDS(sch9_peaks_all, "sch9_peaks_all.rds")

plot <- ggplot(sc_fishers[WT_odds < Inf & WT_odds > 0]) + 
  #geom_rect(aes(xmin=0, xmax=Inf, ymin=-log(0.05,10), ymax=Inf), fill = "#DEEBF7", alpha = 0.3) +
  geom_point(aes(log2(WT_odds), -log10(WT_padj)), color = "gray75", fill = "gray75", alpha = 0.5, size = 2) +
  geom_point(data = WT_peaks_all, aes(log2(WT_odds), -log10(WT_padj), color = "#1F78B4"), alpha = 0.5, size = 2) +
  geom_point(data = WT_peaks_all[WT_D4_pause > WT_D2_pause], aes(log2(WT_odds), -log10(WT_padj), color = "red"), alpha = 0.5, size = 2) +
  scale_x_continuous() + scale_y_continuous(limits = c(0,200)) + 
  scale_color_manual(labels = c("Enhanced","Age-dep"), values = c("#1F78B4", "red"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D12/D1)") +
  theme(legend.position = "none", legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), legend.title=element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_yeast.tiff", plot, width = 6, height = 4, dpi = 300)


### Metagene around stalls ###
WT_stalling_peaks <- readRDS("WT_stalling_peaks.rds")
colnames(WT_stalling_peaks)[colnames(WT_stalling_peaks)=="ID"] <- "peak_ID"
WT_stalling_peaks_ageDep <- WT_stalling_peaks[WT_D4_pause > WT_D2_pause]
#saveRDS(WT_stalling_peaks_ageDep, "WT_stalling_peaks_ageDep.rds")
write.csv(WT_stalling_peaks_ageDep, "WT_stalling_peaks_ageDep.csv")
WT_stalling_peaks_dt <- sc_dtA[WT_stalling_peaks[, c(1,6,24:29,48:50)], allow.cartesian = TRUE]
WT_stalling_peaks_dt[, adjusted := position - peak]

test <- WT_stalling_peaks_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("peak_ID"))
setkeyv(WT_stalling_peaks_dt, c("peak_ID"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = peak_ID]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = peak_ID]
test[, WT_D2_1A_rpc_adjusted := mean(WT_D2_1A), by = peak_ID]
test[, WT_D2_2A_rpc_adjusted := mean(WT_D2_2A), by = peak_ID]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = peak_ID]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = peak_ID]
test <- test[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(WT_stalling_peaks_dt$peak_ID, test$peak_ID))
WT_stalling_peaks_dt <- cbind(WT_stalling_peaks_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
WT_stalling_peaks_dt <- cbind(WT_stalling_peaks_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
WT_stalling_peaks_dt <- cbind(WT_stalling_peaks_dt, WT_D2_1A_rpc_adjusted = test[i]$WT_D2_1A_rpc_adjusted)
WT_stalling_peaks_dt <- cbind(WT_stalling_peaks_dt, WT_D2_2A_rpc_adjusted = test[i]$WT_D2_2A_rpc_adjusted)
WT_stalling_peaks_dt <- cbind(WT_stalling_peaks_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
WT_stalling_peaks_dt <- cbind(WT_stalling_peaks_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
WT_stalling_peaks_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
WT_stalling_peaks_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
WT_stalling_peaks_dt[, WT_D2_1A_norm := WT_D2_1A / WT_D2_1A_rpc_adjusted]
WT_stalling_peaks_dt[, WT_D2_2A_norm := WT_D2_2A / WT_D2_2A_rpc_adjusted]
WT_stalling_peaks_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
WT_stalling_peaks_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
WT_stalling_peaks_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
WT_stalling_peaks_dt[, WT_D2_norm := (WT_D2_1A_norm + WT_D2_2A_norm) / 2]
WT_stalling_peaks_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]
WT_stalling_peaks_dt[, position_norm := position / length]
# saveRDS(WT_stalling_peaks_dt, "WT_stalling_peaks_dt.rds")

# Plot
WT_stalling_peaks_dt <- readRDS("WT_stalling_peaks_dt.rds")
plot <- ggplot(data = WT_stalling_peaks_dt[WT_D0_1A_rpc_adjusted >= 1 & WT_D0_2A_rpc_adjusted >= 1 &
                                             WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                          WT_D4_1A_rpc_adjusted >= 1 & WT_D4_2A_rpc_adjusted >= 1 &
                          peak > 25 & (peak - length) < -25,
                        .(adjusted, WT0 = movingAverage(WT_D0_norm, n=1, center=T),
                          WT2 = movingAverage(WT_D2_norm, n=1, center=T),
                          WT4 = movingAverage(WT_D4_norm, n=1, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'black') +  
  stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = 'blue') +  
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'red')
  #stat_summary(aes(adjusted, WT0), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
   #     fun.args=list(conf.int=0.5), fill = 'black') +
  #stat_summary(aes(adjusted, WT4), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
   #            fun.args=list(conf.int=0.5), fill = 'red')
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

ggplot(data = WT_stalling_peaks_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                     WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 &
                                     WT_D2_1A_rpc_adjusted >= 0.1 & WT_D2_2A_rpc_adjusted >= 0.1 &
                                     peak > 25 & (peak - length) < -25,
                                   .(adjusted, WT0 = movingAverage(WT_D0_norm, n=1, center=T),
                                     WT2 = movingAverage(WT_D2_norm, n=1, center=T),
                                     WT4 = movingAverage(WT_D4_norm, n=1, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'black') +  
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'red')


WT_stalling_peaks_ageDep_orfs <- WT_stalling_peaks_ageDep[, .SD[which.min(peak)], by = orf]
write.csv(WT_stalling_peaks_ageDep_orfs, "WT_stalling_ageDep_Proteins.csv")
WT_stalling_peaks_ageDep_dt <- sc_dtA[WT_stalling_peaks_ageDep[, c(1,6,24:29,48:50)], allow.cartesian = TRUE]
WT_stalling_peaks_ageDep_dt[, adjusted := position - peak]
test <- WT_stalling_peaks_ageDep_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("peak_ID"))
setkeyv(WT_stalling_peaks_ageDep_dt, c("peak_ID"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = peak_ID]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = peak_ID]
test[, WT_D2_1A_rpc_adjusted := mean(WT_D2_1A), by = peak_ID]
test[, WT_D2_2A_rpc_adjusted := mean(WT_D2_2A), by = peak_ID]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = peak_ID]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = peak_ID]
test[, sch9_D0_1A_rpc_adjusted := mean(sch9_D0_1A), by = peak_ID]
test[, sch9_D0_2A_rpc_adjusted := mean(sch9_D0_2A), by = peak_ID]
test[, sch9_D4_1A_rpc_adjusted := mean(sch9_D4_1A), by = peak_ID]
test[, sch9_D4_2A_rpc_adjusted := mean(sch9_D4_2A), by = peak_ID]
test <- test[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(WT_stalling_peaks_ageDep_dt$peak_ID, test$peak_ID))
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, WT_D2_1A_rpc_adjusted = test[i]$WT_D2_1A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, WT_D2_2A_rpc_adjusted = test[i]$WT_D2_2A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, sch9_D0_1A_rpc_adjusted = test[i]$sch9_D0_1A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, sch9_D0_2A_rpc_adjusted = test[i]$sch9_D0_2A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, sch9_D4_1A_rpc_adjusted = test[i]$sch9_D4_1A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt <- cbind(WT_stalling_peaks_ageDep_dt, sch9_D4_2A_rpc_adjusted = test[i]$sch9_D4_2A_rpc_adjusted)
WT_stalling_peaks_ageDep_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, WT_D2_1A_norm := WT_D2_1A / WT_D2_1A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, WT_D2_2A_norm := WT_D2_2A / WT_D2_2A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, sch9_D0_1A_norm := sch9_D0_1A / sch9_D0_1A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, sch9_D0_2A_norm := sch9_D0_2A / sch9_D0_2A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, sch9_D4_1A_norm := sch9_D4_1A / sch9_D4_1A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, sch9_D4_2A_norm := sch9_D4_2A / sch9_D4_2A_rpc_adjusted]
WT_stalling_peaks_ageDep_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
WT_stalling_peaks_ageDep_dt[, WT_D2_norm := (WT_D2_1A_norm + WT_D2_2A_norm) / 2]
WT_stalling_peaks_ageDep_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]
WT_stalling_peaks_ageDep_dt[, sch9_D0_norm := (sch9_D0_1A_norm + sch9_D0_2A_norm) / 2]
WT_stalling_peaks_ageDep_dt[, sch9_D4_norm := (sch9_D4_1A_norm + sch9_D4_2A_norm) / 2]
WT_stalling_peaks_ageDep_dt[, position_norm := position / length]
# saveRDS(WT_stalling_peaks_ageDep_dt, "WT_stalling_peaks_ageDep_dt.rds")

# Plot
WT_stalling_peaks_ageDep_dt <- readRDS("WT_stalling_peaks_ageDep_dt.rds")
ggplot(data = WT_stalling_peaks_ageDep_dt[WT_D0_1A_rpc_adjusted >= 1 & WT_D0_2A_rpc_adjusted >= 1 & 
                                            WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                                            WT_D4_1A_rpc_adjusted >= 1 & WT_D4_2A_rpc_adjusted >= 1 &
                                            peak > 25 & (peak - length) < -25,
                                          .(adjusted, WT0 = movingAverage(WT_D0_norm, n=2, center=T),
                                            WT2 = movingAverage(WT_D2_norm, n=2, center=T),
                                            WT4 = movingAverage(WT_D4_norm, n=2, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = 'blue') +  
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'red') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'black')
#stat_summary(aes(adjusted, WT0), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
#     fun.args=list(conf.int=0.5), fill = 'black') +
#stat_summary(aes(adjusted, WT4), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
#            fun.args=list(conf.int=0.5), fill = 'red')




### Stalling at proline ###
temp <- sc_dtA[motif2 == "PP" & WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 &
                 WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5, c(1:2,6)]
PP_dt <- sc_dtA[temp, allow.cartesian=TRUE]
PP_dt[, adjusted := position - i.position]

# Normalize across proline interval
test <- PP_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("i.ID"))
setkeyv(PP_dt, c("i.ID"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = i.ID]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = i.ID]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = i.ID]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = i.ID]
test <- test[, .SD[which.min(position)], by = i.ID]
i <- cbind(match(PP_dt$i.ID, test$i.ID))
PP_dt <- cbind(PP_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
PP_dt <- cbind(PP_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
PP_dt <- cbind(PP_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
PP_dt <- cbind(PP_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
PP_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
PP_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
PP_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
PP_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
# saveRDS(PP_dt, "PP_dt.rds")


# Plot ribosome occupancy around site of interest
PP_dt <- readRDS("PP_dt.rds")

plot <- ggplot(data = PP_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                      WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 &
                      i.position > 70 & (i.position - length) < -25,
                    .(adjusted, D1 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                      D12 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=1.25, color = 'black') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=1.25, color = 'red')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from xPP site (codons)") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/prolineMetagene.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Gene ontology
GO <- read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/GO_stalling.csv", stringsAsFactors = T, header = T)
GO <- as.data.table(GO)
plot <- ggplot(GO[Organism == "yeast"], aes(x = reorder(Description, Fold.Enrichment), y = Fold.Enrichment, fill = factor(Category))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("BP","CC"), values = c("#A6CEE3", "#B2DF8A"), name = "")
plot <- plot + theme_classic(12) + labs(y = "Fold enrichment", x = "") +
  theme(legend.position = c(0.8,0.2), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_stallingYeast1.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Isolate stall sites using Fisher's exact test ###
sc_dtA <- readRDS("sc_dtA.rds")
sch9_fishers_rep <- sc_dtA[sch9_D0_1A_rpc >= 0.5 & sch9_D0_2A_rpc >= 0.5 &
                             sch9_D0_1A_sum >= 64 & sch9_D0_2A_sum >= 64 &
                             sch9_D2_1A_rpc >= 0.5 & sch9_D2_2A_rpc >= 0.5 &
                             sch9_D2_1A_sum >= 64 & sch9_D2_2A_sum >= 64 &
                             sch9_D4_1A_rpc >= 0.5 & sch9_D4_2A_rpc >= 0.5 &
                             sch9_D4_1A_sum >= 64 & sch9_D4_2A_sum >= 64]
sch9_fishers_rep[, sch9_D0 := round((sch9_D0_1A + sch9_D0_2A) / 2)]
sch9_fishers_rep[, sch9_D2 := round((sch9_D2_1A + sch9_D2_2A) / 2)]
sch9_fishers_rep[, sch9_D4 := round((sch9_D4_1A + sch9_D4_2A) / 2)]
sch9_fishers_rep[, sch9_D0_rpc := round((sch9_D0_1A_rpc + sch9_D0_2A_rpc) / 2)]
sch9_fishers_rep[, sch9_D2_rpc := round((sch9_D2_1A_rpc + sch9_D2_2A_rpc) / 2)]
sch9_fishers_rep[, sch9_D4_rpc := round((sch9_D4_1A_rpc + sch9_D4_2A_rpc) / 2)]
sch9_fishers_rep[, sch9_D0_sum := round((sch9_D0_1A_sum + sch9_D0_2A_sum) / 2)]
sch9_fishers_rep[, sch9_D2_sum := round((sch9_D2_1A_sum + sch9_D2_2A_sum) / 2)]
sch9_fishers_rep[, sch9_D4_sum := round((sch9_D4_1A_sum + sch9_D4_2A_sum) / 2)]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D0_1A_sum - sch9_D0_1A) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D0_2A_sum - sch9_D0_2A) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D2_1A_sum - sch9_D2_1A) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D2_2A_sum - sch9_D2_2A) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D4_1A_sum - sch9_D4_1A) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D4_2A_sum - sch9_D4_2A) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D0_sum - sch9_D0) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D2_sum - sch9_D2) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[!sch9_fishers_rep$orf %in% sch9_fishers_rep[(sch9_D4_sum - sch9_D4) < 0]$orf]
sch9_fishers_rep <- sch9_fishers_rep[, c(1:7,14:19,52,57,62,67,72,77,112:120)]

for (i in 1:nrow(sch9_fishers_rep)) {
  print(i)
  countsD0 <- matrix(c(sch9_fishers_rep[i]$sch9_D0_1A, sch9_fishers_rep[i]$sch9_D0_2A, 
                       (sch9_fishers_rep[i]$sch9_D0_1A_sum - sch9_fishers_rep[i]$sch9_D0_1A), 
                       (sch9_fishers_rep[i]$sch9_D0_2A_sum - sch9_fishers_rep[i]$sch9_D0_2A)), nrow = 2)
  countsD2 <- matrix(c(sch9_fishers_rep[i]$sch9_D2_1A, sch9_fishers_rep[i]$sch9_D2_2A, 
                       (sch9_fishers_rep[i]$sch9_D2_1A_sum - sch9_fishers_rep[i]$sch9_D2_1A), 
                       (sch9_fishers_rep[i]$sch9_D2_2A_sum - sch9_fishers_rep[i]$sch9_D2_2A)), nrow = 2)
  countsD4 <- matrix(c(sch9_fishers_rep[i]$sch9_D4_1A, sch9_fishers_rep[i]$sch9_D4_2A, 
                       (sch9_fishers_rep[i]$sch9_D4_1A_sum - sch9_fishers_rep[i]$sch9_D4_1A), 
                       (sch9_fishers_rep[i]$sch9_D4_2A_sum - sch9_fishers_rep[i]$sch9_D4_2A)), nrow = 2)
  countsD0_D2 <- matrix(c(sch9_fishers_rep[i]$sch9_D2, sch9_fishers_rep[i]$sch9_D0, 
                          (sch9_fishers_rep[i]$sch9_D2_sum - sch9_fishers_rep[i]$sch9_D2), 
                          (sch9_fishers_rep[i]$sch9_D0_sum - sch9_fishers_rep[i]$sch9_D0)), nrow = 2)
  countsD2_D4 <- matrix(c(sch9_fishers_rep[i]$sch9_D4, sch9_fishers_rep[i]$sch9_D2, 
                          (sch9_fishers_rep[i]$sch9_D4_sum - sch9_fishers_rep[i]$sch9_D4), 
                          (sch9_fishers_rep[i]$sch9_D2_sum - sch9_fishers_rep[i]$sch9_D2)), nrow = 2)
  stalling_testD0 <- fisher.test(countsD0)
  sch9_fishers_rep[i, D0_odds := stalling_testD0$estimate]
  sch9_fishers_rep[i, D0_pvalue := stalling_testD0$p.value]
  stalling_testD2 <- fisher.test(countsD2)
  sch9_fishers_rep[i, D2_odds := stalling_testD2$estimate]
  sch9_fishers_rep[i, D2_pvalue := stalling_testD2$p.value]
  stalling_testD4 <- fisher.test(countsD4)
  sch9_fishers_rep[i, D4_odds := stalling_testD4$estimate]
  sch9_fishers_rep[i, D4_pvalue := stalling_testD4$p.value]
  stalling_testD0_D2 <- fisher.test(countsD0_D2)
  sch9_fishers_rep[i, D0_D2_odds := stalling_testD0_D2$estimate]
  sch9_fishers_rep[i, D0_D2_pvalue := stalling_testD0_D2$p.value]
  stalling_testD2_D4 <- fisher.test(countsD2_D4)
  sch9_fishers_rep[i, D4_D2_odds := stalling_testD2_D4$estimate]
  sch9_fishers_rep[i, D4_D2_pvalue := stalling_testD2_D4$p.value]
}
sch9_fishers_rep[, D0_padj := p.adjust(D0_pvalue, method = "BH"), by = orf]
sch9_fishers_rep[, D2_padj := p.adjust(D2_pvalue, method = "BH"), by = orf]
sch9_fishers_rep[, D4_padj := p.adjust(D4_pvalue, method = "BH"), by = orf]
sch9_fishers_rep[, D0_D2_padj := p.adjust(D0_D2_pvalue, method = "BH"), by = orf]
sch9_fishers_rep[, D4_D2_padj := p.adjust(D4_D2_pvalue, method = "BH"), by = orf]
#saveRDS(sc_fishersD0, "sc_fishersD0.rds")


plot <- ggplot(sch9_fishers_rep[D0_odds < Inf & D0_odds > 0]) + 
  geom_point(aes(log2(D0_odds), -log10(D0_padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = sch9_fishers_rep[D0_odds > -Inf & D0_odds < Inf & D0_padj < 0.05 & 
                                       position > 20 & stopdist < -20 &
                                       sch9_D0 > sch9_D0_rpc], aes(log2(D0_odds), -log10(D0_padj)), color = "#E31A1C", alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous(limits = c(0,200)) + 
  scale_color_manual(labels = c("Differential pausing"), values = c("#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D0 Rep1 / D0 Rep2)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_sch9_D0.tiff", plot, width = 6, height = 4, dpi = 300)

plot <- ggplot(sch9_fishers_rep[D2_odds < Inf & D2_odds > 0]) + 
  geom_point(aes(log2(D2_odds), -log10(D2_padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = sch9_fishers_rep[D2_odds > -Inf & D2_odds < Inf & D2_padj < 0.05 & 
                                       position > 20 & stopdist < -20 &
                                       sch9_D2 > sch9_D2_rpc], aes(log2(D2_odds), -log10(D2_padj)), color = "#E31A1C", alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous(limits = c(0,200)) + 
  scale_color_manual(labels = c("Differential pausing"), values = c("#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D2 Rep1 / D2 Rep2)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_sch9_D2.tiff", plot, width = 6, height = 4, dpi = 300)

plot <- ggplot(sch9_fishers_rep[D4_odds < Inf & D4_odds > 0]) + 
  geom_point(aes(log2(D4_odds), -log10(D4_padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = sch9_fishers_rep[D4_odds > -Inf & D4_odds < Inf & D4_padj < 0.05 & 
                                       position > 20 & stopdist < -20 &
                                       sch9_D4 > sch9_D4_rpc], aes(log2(D4_odds), -log10(D4_padj)), color = "#E31A1C", alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous(limits = c(0,200)) + 
  scale_color_manual(labels = c("Differential pausing"), values = c("#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D4 Rep1 / D4 Rep2)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_sch9_D4.tiff", plot, width = 6, height = 4, dpi = 300)

plot <- ggplot(sch9_fishers_rep[D0_D2_odds < Inf & D0_D2_odds > 0]) + 
  geom_point(aes(log2(D0_D2_odds), -log10(D0_D2_padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = sch9_fishers_rep[D0_D2_odds > 1 & D0_D2_odds < Inf & D0_D2_padj < 0.05 & 
                                       position > 20 & stopdist < -20 &
                                       sch9_D2 > sch9_D2_rpc], aes(log2(D0_D2_odds), -log10(D0_D2_padj)), color = "#E31A1C", alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous(limits = c(0,200)) + 
  scale_color_manual(labels = c("Differential pausing"), values = c("#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D2/D0)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_sch9_D0_D2.tiff", plot, width = 6, height = 4, dpi = 300)

plot <- ggplot(sch9_fishers_rep[D4_D2_odds < Inf & D4_D2_odds > 0]) + 
  geom_point(aes(log2(D4_D2_odds), -log10(D4_D2_padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = sch9_fishers_rep[D4_D2_odds > 1 & D4_D2_odds < Inf & D4_D2_padj < 0.05 & 
                                       position > 20 & stopdist < -20 &
                                       sch9_D4 > sch9_D4_rpc], aes(log2(D4_D2_odds), -log10(D4_D2_padj)), color = "#E31A1C", alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous(limits = c(0,200)) + 
  scale_color_manual(labels = c("Differential pausing"), values = c("#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D4/D2)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_sch9_D4_D2.tiff", plot, width = 6, height = 4, dpi = 300)



### Age-independent comparison ###
WT_stalling_peaks <- readRDS("WT_stalling_peaks.rds")
WT_stalling_peaks_ageDep <- readRDS("WT_stalling_peaks_ageDep.rds")
colnames(WT_stalling_peaks)[colnames(WT_stalling_peaks)=="ID"] <- "peak_ID"

WT_stalling_peaks_ageIndep <- WT_stalling_peaks[WT_D4_pause < WT_D2_pause]
WT_stalling_peaks_ageIndep_orfs <- WT_stalling_peaks_ageIndep[, .SD[which.min(peak)], by = orf]
WT_stalling_peaks_ageIndep_dt <- sc_dtA[WT_stalling_peaks_ageIndep[, c(1,6,24:29,48:50)], allow.cartesian = TRUE]
WT_stalling_peaks_ageIndep_dt[, adjusted := position - peak]
test <- WT_stalling_peaks_ageIndep_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("peak_ID"))
setkeyv(WT_stalling_peaks_ageIndep_dt, c("peak_ID"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = peak_ID]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = peak_ID]
test[, WT_D2_1A_rpc_adjusted := mean(WT_D2_1A), by = peak_ID]
test[, WT_D2_2A_rpc_adjusted := mean(WT_D2_2A), by = peak_ID]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = peak_ID]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = peak_ID]
test <- test[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(WT_stalling_peaks_ageIndep_dt$peak_ID, test$peak_ID))
WT_stalling_peaks_ageIndep_dt <- cbind(WT_stalling_peaks_ageIndep_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
WT_stalling_peaks_ageIndep_dt <- cbind(WT_stalling_peaks_ageIndep_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
WT_stalling_peaks_ageIndep_dt <- cbind(WT_stalling_peaks_ageIndep_dt, WT_D2_1A_rpc_adjusted = test[i]$WT_D2_1A_rpc_adjusted)
WT_stalling_peaks_ageIndep_dt <- cbind(WT_stalling_peaks_ageIndep_dt, WT_D2_2A_rpc_adjusted = test[i]$WT_D2_2A_rpc_adjusted)
WT_stalling_peaks_ageIndep_dt <- cbind(WT_stalling_peaks_ageIndep_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
WT_stalling_peaks_ageIndep_dt <- cbind(WT_stalling_peaks_ageIndep_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
WT_stalling_peaks_ageIndep_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
WT_stalling_peaks_ageIndep_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
WT_stalling_peaks_ageIndep_dt[, WT_D2_1A_norm := WT_D2_1A / WT_D2_1A_rpc_adjusted]
WT_stalling_peaks_ageIndep_dt[, WT_D2_2A_norm := WT_D2_2A / WT_D2_2A_rpc_adjusted]
WT_stalling_peaks_ageIndep_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
WT_stalling_peaks_ageIndep_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
WT_stalling_peaks_ageIndep_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
WT_stalling_peaks_ageIndep_dt[, WT_D2_norm := (WT_D2_1A_norm + WT_D2_2A_norm) / 2]
WT_stalling_peaks_ageIndep_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]
WT_stalling_peaks_ageIndep_dt[, position_norm := position / length]
# saveRDS(WT_stalling_peaks_ageIndep_dt, "WT_stalling_peaks_ageIndep_dt.rds")

# Plot
WT_stalling_peaks_dt <- readRDS("WT_stalling_peaks_dt.rds")
WT_stalling_peaks_ageDep_dt <- readRDS("WT_stalling_peaks_ageDep_dt.rds")
plot <- ggplot(data = WT_stalling_peaks_ageIndep_dt[WT_D0_1A_rpc_adjusted >= 1 & WT_D0_2A_rpc_adjusted >= 1 & 
                                                      WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                                                      WT_D4_1A_rpc_adjusted >= 1 & WT_D4_2A_rpc_adjusted >= 1 &
                                                      peak > 25 & (peak - length) < -25,
                                                    .(adjusted, WT0 = movingAverage(WT_D0_norm, n=2, center=T),
                                                      WT2 = movingAverage(WT_D2_norm, n=2, center=T),
                                                      WT4 = movingAverage(WT_D4_norm, n=2, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = '#FB9A99') +  
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E31A1C') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray30') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray30", "#FB9A99", "#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_yeast_ageIndep.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = WT_stalling_peaks_ageDep_dt[WT_D0_1A_rpc_adjusted >= 1 & WT_D0_2A_rpc_adjusted >= 1 & 
                                                    WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                                                    WT_D4_1A_rpc_adjusted >= 1 & WT_D4_2A_rpc_adjusted >= 1 &
                                                    peak > 25 & (peak - length) < -25,
                                                  .(adjusted, WT0 = movingAverage(WT_D0_norm, n=2, center=T),
                                                    WT2 = movingAverage(WT_D2_norm, n=2, center=T),
                                                    WT4 = movingAverage(WT_D4_norm, n=2, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = '#FB9A99') +  
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E31A1C') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray30') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray30", "#FB9A99", "#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_yeast_ageDep.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = WT_stalling_peaks_dt[WT_D0_1A_rpc_adjusted >= 1 & WT_D0_2A_rpc_adjusted >= 1 & 
                                             WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                                             WT_D4_1A_rpc_adjusted >= 1 & WT_D4_2A_rpc_adjusted >= 1 &
                                             peak > 25 & (peak - length) < -25,
                                           .(adjusted, WT0 = movingAverage(WT_D0_norm, n=2, center=T),
                                             WT2 = movingAverage(WT_D2_norm, n=2, center=T),
                                             WT4 = movingAverage(WT_D4_norm, n=2, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = '#FB9A99') +  
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E31A1C') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray30') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray30", "#FB9A99", "#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_yeast_all.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Cys and Trp under 100 in A-site
logomaker(WT_stalling_peaks_ageDep[WT_D4_pause > 6 & residue != "C" & residue != "W"]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$motif3, type = "Logo", bg = bg, color_seed = 6)

logomaker(WT_stalling_peaks_ageDep[WT_D4_pause > 6 & residue != "C" & residue != "W"]$motif3, type = "EDLogo", bg = bg, color_seed = 6)

logomaker(WT_stalling_peaks_ageIndep[WT_D4_pause > 6]$motif3, type = "Logo", bg = bg, color_seed = 6)
logomaker(WT_stalling_peaks[WT_D4_pause > 6]$motif3, type = "Logo", bg = bg, color_seed = 6)
logomaker(WT_stalling_peaks[WT_D2_pause > 6]$motif3, type = "EDLogo", bg = bg, color_seed = 6)



### Faster translation with age
WT_fast_all <- sc_fishers[WT_odds < 1 & WT_odds > -Inf & WT_padj < 0.05 & 
                            position > 20 & stopdist < -20 &
                            WT_D0 > WT_D0_rpc]

orfs <- WT_fast_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- WT_fast_all[g]$WT_odds
  peaks <- WT_fast_all[g]$position
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
WT_D0_peaks <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                          peak = final_peaks,
                          WT_D4_odds = final_odds)
WT_D0_peaks[, ID := as.character(base::paste(orf, peak, sep = "_"))]
setkeyv(WT_D0_peaks, c("orf"))
WT_fast_peaks <- WT_fast_all[WT_fast_all$ID %in% WT_D0_peaks$ID]
i <- cbind(match(WT_fast_peaks$ID, WT_D0_peaks$ID))
WT_fast_peaks <- cbind(WT_fast_peaks, peak = WT_D0_peaks[i]$peak)
WT_fast_peaks <- cbind(WT_fast_peaks, peak_start = WT_D0_peaks[i]$peak_start)
WT_fast_peaks <- cbind(WT_fast_peaks, peak_end = WT_D0_peaks[i]$peak_end)
i <- cbind(match(WT_fast_peaks$ID, sc_dtA$ID))
WT_fast_peaks <- cbind(WT_fast_peaks, motif3 = sc_dtA[i]$motif3)
WT_fast_peaks <- cbind(WT_fast_peaks, motif2 = sc_dtA[i]$motif2)

colnames(WT_fast_peaks)[colnames(WT_fast_peaks)=="ID"] <- "peak_ID"
WT_fast_peaks_dt <- sc_dtA[WT_fast_peaks[, c(1,6,24:29,48:50)], allow.cartesian = TRUE]
WT_fast_peaks_dt[, adjusted := position - peak]

test <- WT_fast_peaks_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("peak_ID"))
setkeyv(WT_fast_peaks_dt, c("peak_ID"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = peak_ID]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = peak_ID]
test[, WT_D2_1A_rpc_adjusted := mean(WT_D2_1A), by = peak_ID]
test[, WT_D2_2A_rpc_adjusted := mean(WT_D2_2A), by = peak_ID]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = peak_ID]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = peak_ID]
test <- test[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(WT_fast_peaks_dt$peak_ID, test$peak_ID))
WT_fast_peaks_dt <- cbind(WT_fast_peaks_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
WT_fast_peaks_dt <- cbind(WT_fast_peaks_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
WT_fast_peaks_dt <- cbind(WT_fast_peaks_dt, WT_D2_1A_rpc_adjusted = test[i]$WT_D2_1A_rpc_adjusted)
WT_fast_peaks_dt <- cbind(WT_fast_peaks_dt, WT_D2_2A_rpc_adjusted = test[i]$WT_D2_2A_rpc_adjusted)
WT_fast_peaks_dt <- cbind(WT_fast_peaks_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
WT_fast_peaks_dt <- cbind(WT_fast_peaks_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
WT_fast_peaks_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
WT_fast_peaks_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
WT_fast_peaks_dt[, WT_D2_1A_norm := WT_D2_1A / WT_D2_1A_rpc_adjusted]
WT_fast_peaks_dt[, WT_D2_2A_norm := WT_D2_2A / WT_D2_2A_rpc_adjusted]
WT_fast_peaks_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
WT_fast_peaks_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
WT_fast_peaks_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
WT_fast_peaks_dt[, WT_D2_norm := (WT_D2_1A_norm + WT_D2_2A_norm) / 2]
WT_fast_peaks_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]
WT_fast_peaks_dt[, position_norm := position / length]

plot <- ggplot(data = WT_fast_peaks_dt[WT_D0_1A_rpc_adjusted >= 1 & WT_D0_2A_rpc_adjusted >= 1 & 
                                         WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                                         WT_D4_1A_rpc_adjusted >= 1 & WT_D4_2A_rpc_adjusted >= 1 &
                                         peak > 25 & (peak - length) < -25,
                                       .(adjusted, WT0 = movingAverage(WT_D0_norm, n=2, center=T),
                                         WT2 = movingAverage(WT_D2_norm, n=2, center=T),
                                         WT4 = movingAverage(WT_D4_norm, n=2, center=T))]) +
  stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = '#F2CF8C') +  
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray40", "#F2CF8C", "#E7A427"), name = "") +
  scale_x_continuous(expand = expand_scale(), limits = c(-25,25))
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from fast site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/fastertranslationMetagene_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

logomaker(WT_fast_peaks$motif3, type = "EDLogo", bg = bg, color_seed = 6)
