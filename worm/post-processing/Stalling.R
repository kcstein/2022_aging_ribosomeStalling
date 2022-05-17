library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
N2_dtA <- readRDS("N2_dtA.rds")
stalling_peaks <- readRDS("stalling_peaks.rds")
N2_fishers <- readRDS("N2_fishers.rds")


### Tests for overall pausing ###
ggplot(N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & position > 5 & stopdist < -5 & residue != "X"], aes(residue, D1_pause)) + geom_violin(aes(residue, D12_pause), color = "red", fill = NA) + geom_violin(fill = NA) +
  scale_y_log10() + coord_cartesian(ylim = c(0.01,400))
ggplot(N2_dtA[D1_2A_rpc >= 1 & D12_2A_rpc >= 1 & residue != "X"], aes(residue, D1_2A_pause)) + geom_violin(aes(residue, D12_2A_pause), color = "red", fill = NA) + geom_violin(fill = NA) +
  scale_y_log10() + coord_cartesian(ylim = c(0.01,400))
ggplot(N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & position > 5 & stopdist < -5 & residue != "X" &
                D1_1A > 5 & D1_2A > 5 & D12_1A > 5 & D12_2A > 5], aes(residue, (D12_pause / D1_pause))) + geom_violin(fill = NA) +
  scale_y_log10() + coord_cartesian(ylim = c(0.01,300))

ggplot(N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & position > 5 & stopdist < -5 & residue == "R" &
                D1_1A > 5 & D1_2A > 5 & D12_1A > 5 & D12_2A > 5], aes("D1_1A", D1_1A_pause)) + geom_boxplot() +
  geom_boxplot(aes("D12_1A", D12_1A_pause)) + geom_boxplot(aes("D12_2A", D12_2A_pause)) + 
  geom_boxplot(aes("D1_2A", D1_2A_pause)) + scale_y_log10()
summary(N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & position > 5 & stopdist < -5 & residue == "R" &
                 D1_1A > 5 & D1_2A > 5 & D12_1A > 5 & D12_2A > 5]$D12_2A_pause)
wilcox.test(N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & position > 5 & stopdist < -5 & residue == "R" &
                     D1_1A > 2 & D1_2A > 2 & D12_1A > 2 & D12_2A > 2]$D1_2A_pause,
            N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & position > 5 & stopdist < -5 & residue == "R" &
                     D1_1A > 2 & D1_2A > 2 & D12_1A > 2 & D12_2A > 2]$D12_2A_pause)

ggplot(N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & 
                D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                N2_D1_2A_rpc >= 0.5 & N2_D12_2A_rpc >= 0.5 & N2_D6_2A_rpc >= 0.5 &
                D1_1A_sum >= 64 & D1_2A_sum >= 64 & N2_D1_2A_sum >= 64 &
                D12_1A_sum >= 64 & D12_2A_sum >= 64 & N2_D12_2A_sum >= 64 &
                N2_D6_2A_sum >= 64 &
                position > 100 & stopdist < -30], 
       aes(D1_1A_pause)) + stat_ecdf(geom = "step", color = "black") +
  scale_x_log10(limits = c(0.05,30)) +
  stat_ecdf(aes(D1_2A_pause), geom = "step", color = "black") +
  stat_ecdf(aes(D12_1A_pause), geom = "step", color = "red") +
  stat_ecdf(aes(D12_2A_pause), geom = "step", color = "red") +
  stat_ecdf(aes(N2_D12_2A_pause), geom = "step", color = "orange") +
  stat_ecdf(aes(N2_D1_2A_pause), geom = "step", color = "gray") +
  stat_ecdf(aes(N2_D6_2A_pause), geom = "step", color = "blue")

ggplot(N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & 
                D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                N2_D1_2A_rpc >= 0.5 & N2_D12_2A_rpc >= 0.5 & N2_D6_2A_rpc >= 0.5 &
                D1_1A_sum >= 64 & D1_2A_sum >= 64 & N2_D1_2A_sum >= 64 &
                D12_1A_sum >= 64 & D12_2A_sum >= 64 & N2_D12_2A_sum >= 64 &
                N2_D6_2A_sum >= 64 & D1_1A_coverage > 0.7 & D1_2A_coverage > 0.7 &
                D12_1A_coverage > 0.7 & D12_2A_coverage > 0.7 &
                N2_D6_2A_coverage > 0.7 & 
                position > 100 & stopdist < -30], 
       aes(D1_pause)) + stat_ecdf(geom = "step", color = "black") +
  scale_x_log10(limits = c(0.05,30)) +
  stat_ecdf(aes(N2_D6_2A_pause), geom = "step", color = "blue") +
  stat_ecdf(aes(D12_pause), geom = "step", color = "red")

ggplot(N2_dtA[D1_1A_rpc >= 0.5 & D12_1A_rpc >= 0.5 & 
                D1_2A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                N2_D6_2A_rpc >= 0.5 & N2_D6_2A_sum >= 64 &
                D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
                position > 100 & stopdist < -30 & residue == "P"], 
       aes(D1_pause)) + stat_ecdf(geom = "step", color = "black") +
  scale_x_log10(limits = c(0.05,30)) +
  stat_ecdf(aes(N2_D6_2A_pause), geom = "step", color = "blue") +
  stat_ecdf(aes(D12_pause), geom = "step", color = "red")

plot <- ggplot(N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 &
                D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 &
                D1_1A_sum >= 64 & D1_2A_sum >= 64 &
                D12_1A_sum >= 64 & D12_2A_sum >= 64 & position > 100 & stopdist < -30], 
       aes(D1_pause)) + stat_ecdf(geom = "step", color = "black", size = 1.25) +
  stat_ecdf(aes(D12_pause), geom = "step", color = "red", size = 1.25) +
  scale_x_log10(limits = c(0.05,30))
plot <- plot + theme_classic(16) + labs(y = "Cumulative fraction", x = "Pause score") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/globalPausing_worm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Isolate stall sites using Fisher's exact test ###
N2_fishers <- N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 &
                       D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 &
                       D1_1A_sum >= 64 & D1_2A_sum >= 64 &
                       D12_1A_sum >= 64 & D12_2A_sum >= 64]
N2_fishers[, D1 := round((D1_1A + D1_2A) / 2)]
N2_fishers[, D12 := round((D12_1A + D12_2A) / 2)]
N2_fishers[, D1_rpc := round((D1_1A_rpc + D1_2A_rpc) / 2)]
N2_fishers[, D12_rpc := round((D12_1A_rpc + D12_2A_rpc) / 2)]
N2_fishers[, D1_sum := round((D1_1A_sum + D1_2A_sum) / 2)]
N2_fishers[, D12_sum := round((D12_1A_sum + D12_2A_sum) / 2)]
N2_fishers <- N2_fishers[, c(1:11,18:29,78:87)]
N2_fishers <- N2_fishers[!N2_fishers$orf %in% N2_fishers[(D1_sum - D1) < 0]$orf] # Find orfs with position with aberrantly high number of reads
N2_fishers <- N2_fishers[!N2_fishers$orf %in% N2_fishers[(D12_sum - D12) < 0]$orf]
#which(N2_fishers$ID %in% c("R09A8.3.1_825"))
N2_fishers <- N2_fishers[!N2_fishers$orf %in% N2_fishers[(D1_1A_sum - D1_1A) < 0]$orf]
N2_fishers <- N2_fishers[!N2_fishers$orf %in% N2_fishers[(D1_2A_sum - D1_2A) < 0]$orf]
N2_fishers <- N2_fishers[!N2_fishers$orf %in% N2_fishers[(D12_1A_sum - D12_1A) < 0]$orf]
N2_fishers <- N2_fishers[!N2_fishers$orf %in% N2_fishers[(D12_2A_sum - D12_2A) < 0]$orf]
View(N2_fishers[(D1_sum - D1) < 0])
View(N2_fishers[(D12_sum - D12) < 0])

for (i in 1:nrow(N2_fishers)) {
  print(i)
  counts1 <- matrix(c(N2_fishers[i]$D12, N2_fishers[i]$D1, 
                      (N2_fishers[i]$D12_sum - N2_fishers[i]$D12), 
                      (N2_fishers[i]$D1_sum - N2_fishers[i]$D1)), nrow = 2)
  counts2 <- matrix(c(N2_fishers[i]$D12_1A, N2_fishers[i]$D1_1A, 
                      (N2_fishers[i]$D12_1A_sum - N2_fishers[i]$D12_1A), 
                      (N2_fishers[i]$D1_1A_sum - N2_fishers[i]$D1_1A)), nrow = 2)
  counts3 <- matrix(c(N2_fishers[i]$D12_2A, N2_fishers[i]$D1_2A, 
                      (N2_fishers[i]$D12_2A_sum - N2_fishers[i]$D12_2A), 
                      (N2_fishers[i]$D1_2A_sum - N2_fishers[i]$D1_2A)), nrow = 2)
  stalling_test1 <- fisher.test(counts1)
  stalling_test2 <- fisher.test(counts2)
  stalling_test3 <- fisher.test(counts3)
  N2_fishers[i, odds := stalling_test1$estimate]
  N2_fishers[i, pvalue := stalling_test1$p.value]
  N2_fishers[i, odds1 := stalling_test2$estimate]
  N2_fishers[i, pvalue1 := stalling_test2$p.value]
  N2_fishers[i, odds2 := stalling_test3$estimate]
  N2_fishers[i, pvalue2 := stalling_test3$p.value]
}
N2_fishers[, padj := p.adjust(pvalue, method = "BH"), by = orf]
N2_fishers[, padj1 := p.adjust(pvalue1, method = "BH"), by = orf]
N2_fishers[, padj2 := p.adjust(pvalue2, method = "BH"), by = orf]
i <- cbind(match(N2_fishers$ID, N2_dtA$ID))
N2_fishers <- cbind(N2_fishers, N2_D6_2A = N2_dtA[i]$N2_D6_2A)
N2_fishers <- cbind(N2_fishers, N2_D6_2A_sum = N2_dtA[i]$N2_D6_2A_sum)
N2_fishers <- cbind(N2_fishers, N2_D6_2A_rpc = N2_dtA[i]$N2_D6_2A_rpc)
N2_fishers <- cbind(N2_fishers, N2_D6_2A_pause = N2_dtA[i]$N2_D1_2A_pause)
#saveRDS(N2_fishers, "N2_fishers.rds")

stalling_peaks_all <- N2_fishers[odds > 1 & odds < Inf & padj < 0.05 & 
                                   position > 20 & stopdist < -20 &
                                   D12 > D12_rpc]
#saveRDS(stalling_peaks_all, "stalling_peaks_all.rds")
# stalling_peaks_all_replicates <- N2_fishers[odds1 > 1 & odds1 < Inf & padj1 < 0.05 & 
#                                               odds2 > 1 & odds2 < Inf & padj2 < 0.05 & 
#                                               position > 30 & stopdist < -5 &
#                                               D12_1A > D12_1A_rpc & D12_2A > D12_2A_rpc]

plot <- ggplot(N2_fishers[odds < Inf & odds > 0]) + 
  #geom_rect(aes(xmin=0, xmax=Inf, ymin=-log(0.05,10), ymax=Inf), fill = "#DEEBF7", alpha = 0.3) +
  geom_point(aes(log2(odds), -log10(padj)), color = "gray75", fill = "gray75", alpha = 0.5, size = 2) +
  geom_point(data = stalling_peaks_all, aes(log2(odds), -log10(padj), color = "#1F78B4"), fill = "#1F78B4", alpha = 0.8, size = 2) +
  geom_point(data = stalling_peaks_all[D12_pause > N2_D6_2A_pause], aes(log2(odds), -log10(padj), color = "red"), fill = "#1F78B4", alpha = 0.8, size = 2) +
  scale_x_continuous() + scale_y_continuous() + 
  scale_color_manual(labels = c("Enhanced","Age-dep"), values = c("#1F78B4", "red"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D12/D1)") +
  theme(legend.position = "none", legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), legend.title=element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_worms.tiff", plot, width = 6, height = 4, dpi = 300)


### Metagene around stalls ###
colnames(stalling_peaks)[colnames(stalling_peaks)=="ID"] <- "peak_ID"
stalling_peaks_dt <- N2_dtA[stalling_peaks[, c(1,6,24:45)], allow.cartesian = TRUE]
stalling_peaks_dt[, adjusted := position - peak]
test <- stalling_peaks_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("peak_ID"))
setkeyv(stalling_peaks_dt, c("peak_ID"))
test[, D1_1A_rpc_adjusted := mean(D1_1A), by = peak_ID]
test[, D1_2A_rpc_adjusted := mean(D1_2A), by = peak_ID]
test[, D6_2A_rpc_adjusted := mean(N2_D6_2A), by = peak_ID]
test[, D12_1A_rpc_adjusted := mean(D12_1A), by = peak_ID]
test[, D12_2A_rpc_adjusted := mean(D12_2A), by = peak_ID]
test <- test[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(stalling_peaks_dt$peak_ID, test$peak_ID))
stalling_peaks_dt <- cbind(stalling_peaks_dt, D1_1A_rpc_adjusted = test[i]$D1_1A_rpc_adjusted)
stalling_peaks_dt <- cbind(stalling_peaks_dt, D1_2A_rpc_adjusted = test[i]$D1_2A_rpc_adjusted)
stalling_peaks_dt <- cbind(stalling_peaks_dt, D6_2A_rpc_adjusted = test[i]$D6_2A_rpc_adjusted)
stalling_peaks_dt <- cbind(stalling_peaks_dt, D12_1A_rpc_adjusted = test[i]$D12_1A_rpc_adjusted)
stalling_peaks_dt <- cbind(stalling_peaks_dt, D12_2A_rpc_adjusted = test[i]$D12_2A_rpc_adjusted)
stalling_peaks_dt[, D1_1A_norm := D1_1A / D1_1A_rpc_adjusted]
stalling_peaks_dt[, D1_2A_norm := D1_2A / D1_2A_rpc_adjusted]
stalling_peaks_dt[, D6_2A_norm := N2_D6_2A / D6_2A_rpc_adjusted]
stalling_peaks_dt[, D12_1A_norm := D12_1A / D12_1A_rpc_adjusted]
stalling_peaks_dt[, D12_2A_norm := D12_2A / D12_2A_rpc_adjusted]
stalling_peaks_dt[, D1_norm := (D1_1A_norm + D1_2A_norm) / 2]
stalling_peaks_dt[, D12_norm := (D12_1A_norm + D12_2A_norm) / 2]
stalling_peaks_dt[, position_norm := position / length]
#saveRDS(stalling_peaks_dt, "stalling_peaks_dt.rds")

stalling_peaks_dt <- readRDS("stalling_peaks_dt.rds")
plot <- ggplot(stalling_peaks_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                                   D6_2A_rpc_adjusted >= 1 & 
                               D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 &
                               peak > 25 & (peak - length) < -25,
                             .(adjusted, D1 = movingAverage(D1_norm, n=1, center=T),
                               D6 = movingAverage(D6_2A_norm, n=1, center=T),
                               D12 = movingAverage(D12_norm, n=1, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=1.25, color = 'black') +  
  #stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
  #            fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(adjusted, D6), fun.y = "mean", geom = "line", size=1.25, color = 'blue') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=1.25, color = 'red')
#stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
#            fun.args=list(conf.int=0.5), fill = 'red')
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


stalling_peaks_ageDep <- stalling_peaks[D12_pause > N2_D6_2A_pause]
genenames <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/geneNames.csv', header = TRUE, stringsAsFactors = TRUE))
i <- cbind(match(stalling_peaks_ageDep$orf, genenames$Transcript.stable.ID))
stalling_peaks_ageDep <- cbind(stalling_peaks_ageDep, WBgene = genenames[i]$Gene.stable.ID)
hartl_proteome <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/proteinLists/hartl_proteome.csv', header = TRUE, stringsAsFactors = TRUE))
i <- cbind(match(stalling_peaks_ageDep$WBgene, hartl_proteome$WBgene))
stalling_peaks_ageDep <- cbind(stalling_peaks_ageDep, gene = hartl_proteome[i]$gene)
stalling_peaks_ageDep <- cbind(stalling_peaks_ageDep, uniprot = hartl_proteome[i]$uniprot)
#saveRDS(stalling_peaks_ageDep, "stalling_peaks_ageDep.rds")
write.csv(stalling_peaks_ageDep, "stalling_peaks_ageDep.csv")
stalling_peaks_ageDep_orfs <- stalling_peaks_ageDep[, .SD[which.max(D12_pause)], by = orf]
write.csv(stalling_peaks_ageDep_orfs, "stalling_ageDep_Proteins.csv")


stalling_peaks_ageDep_dt <- N2_dtA[stalling_peaks_ageDep[, c(1,6,24:45)], allow.cartesian = TRUE]
stalling_peaks_ageDep_dt[, adjusted := position - peak]
test <- stalling_peaks_ageDep_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("peak_ID"))
setkeyv(stalling_peaks_ageDep_dt, c("peak_ID"))
test[, D1_1A_rpc_adjusted := mean(D1_1A), by = peak_ID]
test[, D1_2A_rpc_adjusted := mean(D1_2A), by = peak_ID]
test[, D6_2A_rpc_adjusted := mean(N2_D6_2A), by = peak_ID]
test[, D12_1A_rpc_adjusted := mean(D12_1A), by = peak_ID]
test[, D12_2A_rpc_adjusted := mean(D12_2A), by = peak_ID]
test <- test[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(stalling_peaks_ageDep_dt$peak_ID, test$peak_ID))
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D1_1A_rpc_adjusted = test[i]$D1_1A_rpc_adjusted)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D1_2A_rpc_adjusted = test[i]$D1_2A_rpc_adjusted)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D6_2A_rpc_adjusted = test[i]$D6_2A_rpc_adjusted)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D12_1A_rpc_adjusted = test[i]$D12_1A_rpc_adjusted)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D12_2A_rpc_adjusted = test[i]$D12_2A_rpc_adjusted)
stalling_peaks_ageDep_dt[, D1_1A_norm := D1_1A / D1_1A_rpc_adjusted]
stalling_peaks_ageDep_dt[, D1_2A_norm := D1_2A / D1_2A_rpc_adjusted]
stalling_peaks_ageDep_dt[, D6_2A_norm := N2_D6_2A / D6_2A_rpc_adjusted]
stalling_peaks_ageDep_dt[, D12_1A_norm := D12_1A / D12_1A_rpc_adjusted]
stalling_peaks_ageDep_dt[, D12_2A_norm := D12_2A / D12_2A_rpc_adjusted]
stalling_peaks_ageDep_dt[, D1_norm := (D1_1A_norm + D1_2A_norm) / 2]
stalling_peaks_ageDep_dt[, D12_norm := (D12_1A_norm + D12_2A_norm) / 2]
stalling_peaks_ageDep_dt[, position_norm := position / length]
#saveRDS(stalling_peaks_ageDep, "stalling_peaks_ageDep.rds")
#saveRDS(stalling_peaks_ageDep_dt, "stalling_peaks_ageDep_dt.rds")

stalling_peaks_ageDep_dt <- readRDS("stalling_peaks_ageDep_dt.rds")
ggplot(stalling_peaks_ageDep_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                                  D6_2A_rpc_adjusted >= 1 & 
                                  D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 &
                                  peak > 25 & (peak - length) < -25,
                                .(adjusted, D1 = movingAverage(D1_norm, n=2, center=T),
                                  D6 = movingAverage(D6_2A_norm, n=2, center=T),
                                  D12 = movingAverage(D12_norm, n=2, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=1.25, color = 'black') +  
  #stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
  #            fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(adjusted, D6), fun.y = "mean", geom = "line", size=1.25, color = 'blue') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=1.25, color = 'red')
#stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
#            fun.args=list(conf.int=0.5), fill = 'red')


### Peaks by z-scores
D1_Zbackground <- N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & 
                           D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                           position > 20 & stopdist < -20]
D1_Zpeaks <- N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & 
                      D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                      D1_1A_Z > 3.5 & D1_2A_Z > 3.5 & 
                      position > 20 & stopdist < -20]
D6_Zbackground <- N2_dtA[N2_D6_2A_rpc >= 0.5 & 
                           N2_D6_2A_sum >= 64 & 
                           position > 20 & stopdist < -20]
D6_Zpeaks <- N2_dtA[N2_D6_2A_rpc >= 0.5 & 
                      N2_D6_2A_sum >= 64 & 
                      N2_D6_2A_Z > 3.5 & 
                      position > 20 & stopdist < -20]
D12_Zbackground <- N2_dtA[D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                            D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
                            position > 20 & stopdist < -20]
D12_Zpeaks <- N2_dtA[D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                       D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
                       D12_1A_Z > 3.5 & D12_2A_Z > 3.5 & 
                       position > 20 & stopdist < -20]

D1freq <- (summary(D1_Zpeaks$residue) / length(D1_Zpeaks$residue)) / 
  (summary(D1_Zbackground$residue) / length(D1_Zbackground$residue))
D6freq <- (summary(D6_Zpeaks$residue) / length(D6_Zpeaks$residue)) / 
  (summary(D6_Zbackground$residue) / length(D6_Zbackground$residue))
D12freq <- (summary(D12_Zpeaks$residue) / length(D12_Zpeaks$residue)) / 
  (summary(D12_Zbackground$residue) / length(D12_Zbackground$residue))
freq1 <- data.table(age = "D1", residue = names(D1freq), freq = D1freq)
freq6 <- data.table(age = "D6", residue = names(D6freq), freq = D6freq)
freq12 <- data.table(age = "D12", residue = names(D12freq), freq = D12freq)
freq <- rbind(freq1, freq6, freq12)
freq$age <- factor(freq$age, levels = c("D1", "D6", "D12"))

ggplot(freq[residue != "X"], aes(residue,log2(freq), fill = age)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(limits = c("D1", "D6", "D12"), 
                    labels = c("Day 1", "Day 6", "Day 12"),
                    values = c("#1F78B4", "#999999", "#E31A1C"), name = "")

### Ribosome collision positions ###
temp2 <- stalling_peaks_dt[adjusted == -10] # positions 10aa upstream of stall site
temp2a <- temp2[temp2$ID %in% stalling_peaks$peak_ID] # subset upstream positions to those that are also stall sites
collisions2 <- stalling_peaks[stalling_peaks$peak_ID %in% temp2a$peak_ID] # obtain stalls that have upstream stall as well
temp3 <- stalling_peaks_dt[adjusted == -20] # positions 20aa upstream of stall site
temp3a <- temp3[temp3$ID %in% stalling_peaks$peak_ID] # subset 20aa upstream positions to those that are also stall sites
collisions3 <- stalling_peaks[(stalling_peaks$peak_ID %in% temp3a$peak_ID) & (stalling_peaks$peak_ID %in% collisions2$peak_ID)] # obtain stalls that have upstream stall as well

write.csv(collisions2, "ribosomeCollisions2.csv")
write.csv(collisions3, "ribosomeCollisions3.csv")

collisions2_dt <- stalling_peaks_dt[stalling_peaks_dt$peak_ID %in% collisions2$peak_ID]
collisions3_dt <- stalling_peaks_dt[stalling_peaks_dt$peak_ID %in% collisions3$peak_ID]

ggplot(data = collisions3_dt[D1_1A_rpc_adjusted >= 0.5 & D1_2A_rpc_adjusted >= 0.5 & 
                             D12_1A_rpc_adjusted >= 0.5 & D12_2A_rpc_adjusted >= 0.5 & 
                             peak > 25 & (peak - length) < -25,
                           .(adjusted, D1 = movingAverage(D1_norm, n=1, center=T),
                             D12 = movingAverage(D12_norm, n=1, center=T))]) + xlim(-40, 40) +
  stat_summary(aes(adjusted, D1), fun.y = "median", geom = "line", size=0.5, color = 'black') +  
  #stat_summary(aes(adjusted, D1), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
  #            fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(adjusted, D12), fun.y = "median", geom = "line", size=0.5, color = 'red')

library(RWebLogo)
weblogo(collisions2$motif3, open = TRUE, file.out = "collisions2.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
collisions2a <- collisions2[order(D12_pause, decreasing = TRUE),]
weblogo(collisions2a[1:30]$motif3, open = TRUE, file.out = "collisions2.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
weblogo(collisions3$motif3, open = TRUE, file.out = "collisions3.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
collisions3a <- collisions3[order(D12_pause, decreasing = TRUE),]
weblogo(collisions3a[1:25]$motif3, open = TRUE, file.out = "collisions3.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
library(Logolas)
bg <- c(A=0.06325907,C=0.02076066,D=0.05304005,E=0.06526124,F=0.04787855,G=0.05344524,H=0.02281571,I=0.06213785,K=0.06339429,L=0.08643033,M=0.02648252,N=0.04877476,P=0.04892958,Q=0.04081622,R=0.05137781,S=0.08080618,T=0.05891986,V=0.06229459,W=0.01107328,Y=0.03210222)
logomaker(collisions2a[1:100]$motif3, type = "EDLogo", bg = bg)
logomaker(collisions3$motif3, type = "EDLogo", bg = bg)


### Ribosome collision positions ###
stalling_peaks_dt <- readRDS("stalling_peaks_dt.rds")
temp <- stalling_peaks_dt[adjusted >= -15 & adjusted <= -10]
temp1 <- temp[temp$ID %in% stalling_peaks$ID]
collisions <- stalling_peaks_dt[stalling_peaks_dt$peak_ID %in% temp1$peak_ID]
collisions_orfs <- collisions[, .SD[which.min(peak)], by = orf]
i <- cbind(match(collisions_orfs$orf, genenames$Transcript.stable.ID))
collisions_orfs <- cbind(collisions_orfs, WBgene = genenames[i]$Gene.stable.ID)
i <- cbind(match(collisions_orfs$WBgene, hartl_proteome$WBgene))
collisions_orfs <- cbind(collisions_orfs, gene = hartl_proteome[i]$gene)
collisions_orfs <- cbind(collisions_orfs, uniprot = hartl_proteome[i]$uniprot)
write.csv(collisions_orfs, "ribosomeCollisions.csv")

logomaker(collisions1[D12_pause > 6]$motif3, type = "EDLogo", bg = bg, color_seed = 6)

ggplot(data = collisions[D1_1A_rpc_adjusted >= 0.5 & D1_2A_rpc_adjusted >= 0.5 & 
                           D12_1A_rpc_adjusted >= 0.5 & D12_2A_rpc_adjusted >= 0.5 & 
                           peak > 25 & (peak - length) < -25,
                         .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=2, center=T),
                           D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=2, center=T))]) + xlim(-40, 40) +
  #stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=0.5, color = 'red')


### Stalling at proline ###
temp <- N2_dtA[motif2 == "PP" & D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 &
                 D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5, c(1:2,6)]
PP_dt <- N2_dtA[temp, allow.cartesian=TRUE]
PP_dt[, adjusted := position - i.position]

# Normalize across proline interval
test <- PP_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("i.ID"))
setkeyv(PP_dt, c("i.ID"))
test[, D1_1A_rpc_adjusted := mean(D1_1A), by = i.ID]
test[, D1_2A_rpc_adjusted := mean(D1_2A), by = i.ID]
test[, D12_1A_rpc_adjusted := mean(D12_1A), by = i.ID]
test[, D12_2A_rpc_adjusted := mean(D12_2A), by = i.ID]
test <- test[, .SD[which.min(position)], by = i.ID]
i <- cbind(match(PP_dt$i.ID, test$i.ID))
PP_dt <- cbind(PP_dt, D1_1A_rpc_adjusted = test[i]$D1_1A_rpc_adjusted)
PP_dt <- cbind(PP_dt, D1_2A_rpc_adjusted = test[i]$D1_2A_rpc_adjusted)
PP_dt <- cbind(PP_dt, D12_1A_rpc_adjusted = test[i]$D12_1A_rpc_adjusted)
PP_dt <- cbind(PP_dt, D12_2A_rpc_adjusted = test[i]$D12_2A_rpc_adjusted)
PP_dt[, D1_1A_norm := D1_1A / D1_1A_rpc_adjusted]
PP_dt[, D1_2A_norm := D1_2A / D1_2A_rpc_adjusted]
PP_dt[, D12_1A_norm := D12_1A / D12_1A_rpc_adjusted]
PP_dt[, D12_2A_norm := D12_2A / D12_2A_rpc_adjusted]
# saveRDS(PP_dt, "PP_dt.rds")

# Plot ribosome occupancy around site of interest
PP_dt <- readRDS("PP_dt.rds")

ggplot(data = PP_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                      D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 &
                      i.position > 70 & (i.position - length) < -25,
                    .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=1, center=T),
                      D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=1, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=0.5, color = 'red')


### Asymmetry around stalls ###
stalling_peaks_ageDep_dt <- readRDS("stalling_peaks_ageDep_dt.rds")
test_up <- stalling_peaks_ageDep_dt[adjusted < 0 & adjusted > -50]
setkeyv(test_up, c("peak_ID"))
setkeyv(stalling_peaks_ageDep_dt, c("peak_ID"))
test_up[, D1_1A_rpc_up := mean(D1_1A), by = peak_ID]
test_up[, D1_2A_rpc_up := mean(D1_2A), by = peak_ID]
test_up[, D12_1A_rpc_up := mean(D12_1A), by = peak_ID]
test_up[, D12_2A_rpc_up := mean(D12_2A), by = peak_ID]
test_up <- test_up[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(stalling_peaks_ageDep_dt$peak_ID, test_up$peak_ID))
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D1_1A_rpc_up = test_up[i]$D1_1A_rpc_up)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D1_2A_rpc_up = test_up[i]$D1_2A_rpc_up)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D12_1A_rpc_up = test_up[i]$D12_1A_rpc_up)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D12_2A_rpc_up = test_up[i]$D12_2A_rpc_up)
test_down <- stalling_peaks_ageDep_dt[adjusted > 0 & adjusted < 50]
setkeyv(test_down, c("peak_ID"))
setkeyv(stalling_peaks_ageDep_dt, c("peak_ID"))
test_down[, D1_1A_rpc_down := mean(D1_1A), by = peak_ID]
test_down[, D1_2A_rpc_down := mean(D1_2A), by = peak_ID]
test_down[, D12_1A_rpc_down := mean(D12_1A), by = peak_ID]
test_down[, D12_2A_rpc_down := mean(D12_2A), by = peak_ID]
test_down <- test_down[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(stalling_peaks_ageDep_dt$peak_ID, test_down$peak_ID))
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D1_1A_rpc_down = test_down[i]$D1_1A_rpc_down)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D1_2A_rpc_down = test_down[i]$D1_2A_rpc_down)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D12_1A_rpc_down = test_down[i]$D12_1A_rpc_down)
stalling_peaks_ageDep_dt <- cbind(stalling_peaks_ageDep_dt, D12_2A_rpc_down = test_down[i]$D12_2A_rpc_down)

D12_asymmetry <- stalling_peaks_ageDep_dt[adjusted == 0]
D12_asymmetry[, D1_1A_asymmetry := D1_1A_rpc_up / D1_1A_rpc_down]
D12_asymmetry[, D1_2A_asymmetry := D1_2A_rpc_up / D1_2A_rpc_down]
D12_asymmetry[, D12_1A_asymmetry := D12_1A_rpc_up / D12_1A_rpc_down]
D12_asymmetry[, D12_2A_asymmetry := D12_2A_rpc_up / D12_2A_rpc_down]
D12_asymmetry[, D1_asymmetry := (D1_1A_asymmetry + D1_2A_asymmetry) / 2]
D12_asymmetry[, D12_asymmetry := (D12_1A_asymmetry + D12_2A_asymmetry) / 2]
D12_asymmetry[, position_norm := position / length]
setkeyv(D12_asymmetry, c("orf"))

ggplot(D12_asymmetry[D1_1A_rpc_up > 1 & D1_2A_rpc_up > 1 &
                       D12_1A_rpc_up > 1 & D12_2A_rpc_up > 1 &
                       D1_1A_rpc_down > 1 & D1_2A_rpc_down > 1 &
                       D12_1A_rpc_down > 1 & D12_2A_rpc_down > 1 & D12_pause >= 10 &
                       peak > 70 & (peak - length) < -70], aes("D1", D1_asymmetry)) + geom_boxplot() +
  geom_boxplot(aes("D12", D12_asymmetry)) + ylim(0,3)
wilcox.test(D12_asymmetry[D1_1A_rpc_up > 0.5 & D1_2A_rpc_up > 0.5 &
                            D12_1A_rpc_up > 0.5 & D12_2A_rpc_up > 0.5 &
                            D1_1A_rpc_down > 0.5 & D1_2A_rpc_down > 0.5 &
                            D12_1A_rpc_down > 0.5 & D12_2A_rpc_down > 0.5 & D12_pause >= 15 &
                            peak > 70 & (peak - length) < -70]$D1_asymmetry, 
            D12_asymmetry[D1_1A_rpc_up > 0.5 & D1_2A_rpc_up > 0.5 &
                            D12_1A_rpc_up > 0.5 & D12_2A_rpc_up > 0.5 &
                            D1_1A_rpc_down > 0.5 & D1_2A_rpc_down > 0.5 &
                            D12_1A_rpc_down > 0.5 & D12_2A_rpc_down > 0.5 & D12_pause >= 15 &
                            peak > 70 & (peak - length) < -70]$D12_asymmetry, alternative = 't')

wilcox.test(D12_asymmetry[D12_pause >= 10]$D1_asymmetry, 
            D12_asymmetry[D12_pause >= 10]$D12_asymmetry, alternative = 't')

ggplot(D12_asymmetry[D12_pause >= 10 & peak > 70 & (peak - length) < -70], aes("D1", D1_asymmetry)) + geom_boxplot() +
  geom_boxplot(aes("D12", D12_asymmetry)) + ylim(0,3)


### Gene ontology
GO <- read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/GO_stalling.csv", stringsAsFactors = T, header = T)
GO <- as.data.table(GO)
plot <- ggplot(GO[Organism == "worms"], aes(x = reorder(Description, -log10(Benjamini)), y = -log10(Benjamini), fill = factor(Category))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("BP","CC"), values = c("#A6CEE3", "#B2DF8A"), name = "")
plot <- plot + theme_classic(12) + labs(y = "-log10(adjusted p-value)", x = "") +
  theme(legend.position = c(0.8,0.2), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_stallingWorm1.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(GO[Organism == "worms"], aes(x = reorder(Description, Fold.Enrichment), y = Fold.Enrichment, fill = factor(Category))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("BP","CC"), values = c("#A6CEE3", "#B2DF8A"), name = "")
plot <- plot + theme_classic(12) + labs(y = "Fold enrichment", x = "") +
  theme(legend.position = c(0.8,0.2), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_stallingWorm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Volcano of replicates
N2_dtA <- readRDS("N2_dtA.rds")
N2_fishers_rep <- N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 &
                           D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 &
                           D1_1A_sum >= 64 & D1_2A_sum >= 64 &
                           D12_1A_sum >= 64 & D12_2A_sum >= 64]
N2_fishers_rep[, D1 := round((D1_1A + D1_2A) / 2)]
N2_fishers_rep[, D12 := round((D12_1A + D12_2A) / 2)]
N2_fishers_rep[, D1_rpc := round((D1_1A_rpc + D1_2A_rpc) / 2)]
N2_fishers_rep[, D12_rpc := round((D12_1A_rpc + D12_2A_rpc) / 2)]
N2_fishers_rep[, D1_sum := round((D1_1A_sum + D1_2A_sum) / 2)]
N2_fishers_rep[, D12_sum := round((D12_1A_sum + D12_2A_sum) / 2)]
N2_fishers_rep <- N2_fishers_rep[, c(1:11,18:29,78:87)]
N2_fishers_rep <- N2_fishers_rep[!N2_fishers_rep$orf %in% N2_fishers_rep[(D1_sum - D1) < 0]$orf] # Find orfs with position with aberrantly high number of reads
N2_fishers_rep <- N2_fishers_rep[!N2_fishers_rep$orf %in% N2_fishers_rep[(D12_sum - D12) < 0]$orf]
N2_fishers_rep <- N2_fishers_rep[!N2_fishers_rep$orf %in% N2_fishers_rep[(D1_1A_sum - D1_1A) < 0]$orf]
N2_fishers_rep <- N2_fishers_rep[!N2_fishers_rep$orf %in% N2_fishers_rep[(D1_2A_sum - D1_2A) < 0]$orf]
N2_fishers_rep <- N2_fishers_rep[!N2_fishers_rep$orf %in% N2_fishers_rep[(D12_1A_sum - D12_1A) < 0]$orf]
N2_fishers_rep <- N2_fishers_rep[!N2_fishers_rep$orf %in% N2_fishers_rep[(D12_2A_sum - D12_2A) < 0]$orf]
View(N2_fishers_rep[(D1_sum - D1) < 0])
View(N2_fishers_rep[(D12_sum - D12) < 0])

for (i in 1:nrow(N2_fishers_rep)) {
  print(i)
  countsD1 <- matrix(c(N2_fishers_rep[i]$D1_1A, N2_fishers_rep[i]$D1_2A, 
                       (N2_fishers_rep[i]$D1_1A_sum - N2_fishers_rep[i]$D1_1A), 
                       (N2_fishers_rep[i]$D1_2A_sum - N2_fishers_rep[i]$D1_2A)), nrow = 2)
  countsD12 <- matrix(c(N2_fishers_rep[i]$D12_1A, N2_fishers_rep[i]$D12_2A, 
                        (N2_fishers_rep[i]$D12_1A_sum - N2_fishers_rep[i]$D12_1A), 
                        (N2_fishers_rep[i]$D12_2A_sum - N2_fishers_rep[i]$D12_2A)), nrow = 2)stalling_test1 <- fisher.test(counts1)
  stalling_testD1 <- fisher.test(countsD1)
  stalling_testD12 <- fisher.test(countsD12)
  N2_fishers_rep[i, D1_odds := stalling_testD1$estimate]
  N2_fishers_rep[i, D1_pvalue := stalling_testD1$p.value]
  N2_fishers_rep[i, D12_odds := stalling_testD12$estimate]
  N2_fishers_rep[i, D12_pvalue := stalling_testD12$p.value]
}
N2_fishers_rep[, D1_padj := p.adjust(D1_pvalue, method = "BH"), by = orf]
N2_fishers_rep[, D12_padj := p.adjust(D12_pvalue, method = "BH"), by = orf]
#saveRDS(N2_fishers_rep, "N2_fishers_rep.rds")

plot <- ggplot(N2_fishers_rep[D1_odds < Inf & D1_odds > 0]) + 
  geom_point(aes(log2(D1_odds), -log10(D1_padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = N2_fishers_rep[D1_odds > 0 & D1_odds < Inf & D1_padj < 0.05 & 
                                     position > 20 & stopdist < -20 &
                                     D1 > D1_rpc], aes(log2(D1_odds), -log10(D1_padj), color = "#E31A1C"), alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous() + 
  scale_color_manual(labels = c("Differential pausing"), values = c("#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D1 Rep1 / D1 Rep2)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_D1.tiff", plot, width = 6, height = 4, dpi = 300)

plot <- ggplot(N2_fishers_rep[D12_odds < Inf & D12_odds > 0]) + 
  geom_point(aes(log2(D12_odds), -log10(D12_padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = N2_fishers_rep[D12_odds > 0 & D12_odds < Inf & D12_padj < 0.05 & 
                                     position > 20 & stopdist < -20 &
                                     D12 > D12_rpc], aes(log2(D12_odds), -log10(D12_padj), color = "#E31A1C"), alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous() + 
  scale_color_manual(labels = c("Differential pausing"), values = c("#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D12 Rep1 / D12 Rep2)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_D12.tiff", plot, width = 6, height = 4, dpi = 300)


# Metagene of age-dependent pausing sites
stalling_peaks <- readRDS("stalling_peaks.rds")
colnames(stalling_peaks)[colnames(stalling_peaks)=="ID"] <- "peak_ID"
stalling_peaks_ageIndep <- stalling_peaks[D12_pause < N2_D6_2A_pause]
stalling_peaks_ageIndep_dt <- N2_dtA[stalling_peaks_ageIndep[, c(1,6,24:45)], allow.cartesian = TRUE]
stalling_peaks_ageIndep_dt[, adjusted := position - peak]
test <- stalling_peaks_ageIndep_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("peak_ID"))
setkeyv(stalling_peaks_ageIndep_dt, c("peak_ID"))
test[, D1_1A_rpc_adjusted := mean(D1_1A), by = peak_ID]
test[, D1_2A_rpc_adjusted := mean(D1_2A), by = peak_ID]
test[, D6_2A_rpc_adjusted := mean(N2_D6_2A), by = peak_ID]
test[, D12_1A_rpc_adjusted := mean(D12_1A), by = peak_ID]
test[, D12_2A_rpc_adjusted := mean(D12_2A), by = peak_ID]
test <- test[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(stalling_peaks_ageIndep_dt$peak_ID, test$peak_ID))
stalling_peaks_ageIndep_dt <- cbind(stalling_peaks_ageIndep_dt, D1_1A_rpc_adjusted = test[i]$D1_1A_rpc_adjusted)
stalling_peaks_ageIndep_dt <- cbind(stalling_peaks_ageIndep_dt, D1_2A_rpc_adjusted = test[i]$D1_2A_rpc_adjusted)
stalling_peaks_ageIndep_dt <- cbind(stalling_peaks_ageIndep_dt, D6_2A_rpc_adjusted = test[i]$D6_2A_rpc_adjusted)
stalling_peaks_ageIndep_dt <- cbind(stalling_peaks_ageIndep_dt, D12_1A_rpc_adjusted = test[i]$D12_1A_rpc_adjusted)
stalling_peaks_ageIndep_dt <- cbind(stalling_peaks_ageIndep_dt, D12_2A_rpc_adjusted = test[i]$D12_2A_rpc_adjusted)
stalling_peaks_ageIndep_dt[, D1_1A_norm := D1_1A / D1_1A_rpc_adjusted]
stalling_peaks_ageIndep_dt[, D1_2A_norm := D1_2A / D1_2A_rpc_adjusted]
stalling_peaks_ageIndep_dt[, D6_2A_norm := N2_D6_2A / D6_2A_rpc_adjusted]
stalling_peaks_ageIndep_dt[, D12_1A_norm := D12_1A / D12_1A_rpc_adjusted]
stalling_peaks_ageIndep_dt[, D12_2A_norm := D12_2A / D12_2A_rpc_adjusted]
stalling_peaks_ageIndep_dt[, D1_norm := (D1_1A_norm + D1_2A_norm) / 2]
stalling_peaks_ageIndep_dt[, D12_norm := (D12_1A_norm + D12_2A_norm) / 2]
stalling_peaks_ageIndep_dt[, position_norm := position / length]
plot <- ggplot(stalling_peaks_ageIndep_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                                            D6_2A_rpc_adjusted >= 1 & 
                                            D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 &
                                            peak > 25 & (peak - length) < -25,
                                          .(adjusted, D1 = movingAverage(D1_norm, n=2, center=T),
                                            D6 = movingAverage(D6_2A_norm, n=2, center=T),
                                            D12 = movingAverage(D12_norm, n=2, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D6), fun.y = "mean", geom = "line", size=1.25, color = '#FB9A99') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=1.25, color = '#E31A1C') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=1.25, color = 'gray30') +
  scale_color_manual(labels = c("Day 1","Day 6", "Day 12"), values = c("gray30", "#FB9A99", "#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_worms_ageIndep.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

stalling_peaks <- readRDS("stalling_peaks.rds")
plot <- ggplot(stalling_peaks_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                                   D6_2A_rpc_adjusted >= 1 & 
                                   D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 &
                                   peak > 25 & (peak - length) < -25,
                                 .(adjusted, D1 = movingAverage(D1_norm, n=2, center=T),
                                   D6 = movingAverage(D6_2A_norm, n=2, center=T),
                                   D12 = movingAverage(D12_norm, n=2, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D6), fun.y = "mean", geom = "line", size=1.25, color = '#FB9A99') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=1.25, color = '#E31A1C') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=1.25, color = 'gray30') +
  scale_color_manual(labels = c("Day 1","Day 6", "Day 12"), values = c("gray30", "#FB9A99", "#E31A1C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_worms_all.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


stalling_peaks_ageIndep <- readRDS("stalling_peaks_ageIndep.rds")
logomaker(stalling_peaks_ageDep[D12_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(stalling_peaks_ageIndep[D12_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(stalling_peaks[D12_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)



### Fast peaks
fast_peaks_all <- N2_fishers[odds < 1 & odds > -Inf & padj < 0.05 & 
                               position > 20 & stopdist < -20 &
                               D1 > D1_rpc]

orfs <- fast_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- fast_peaks_all[g]$odds
  peaks <- fast_peaks_all[g]$position
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
fast_peaks1 <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                          peak = final_peaks,
                          odds = final_odds)
fast_peaks1[, ID := as.character(base::paste(orf, peak, sep = "_"))]
setkeyv(fast_peaks1, c("orf"))
fast_peaks <- fast_peaks_all[fast_peaks_all$ID %in% fast_peaks1$ID]
i <- cbind(match(fast_peaks$ID, fast_peaks1$ID))
fast_peaks <- cbind(fast_peaks, peak = fast_peaks1[i]$peak)
fast_peaks <- cbind(fast_peaks, peak_start = fast_peaks1[i]$peak_start)
fast_peaks <- cbind(fast_peaks, peak_end = fast_peaks1[i]$peak_end)

colnames(fast_peaks)[colnames(fast_peaks)=="ID"] <- "peak_ID"
fast_peaks_dt <- N2_dtA[fast_peaks[, c(1,6,24:45)], allow.cartesian = TRUE]
fast_peaks_dt[, adjusted := position - peak]
test <- fast_peaks_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("peak_ID"))
setkeyv(fast_peaks_dt, c("peak_ID"))
test[, D1_1A_rpc_adjusted := mean(D1_1A), by = peak_ID]
test[, D1_2A_rpc_adjusted := mean(D1_2A), by = peak_ID]
test[, D6_2A_rpc_adjusted := mean(N2_D6_2A), by = peak_ID]
test[, D12_1A_rpc_adjusted := mean(D12_1A), by = peak_ID]
test[, D12_2A_rpc_adjusted := mean(D12_2A), by = peak_ID]
test <- test[, .SD[which.min(position)], by = peak_ID]
i <- cbind(match(fast_peaks_dt$peak_ID, test$peak_ID))
fast_peaks_dt <- cbind(fast_peaks_dt, D1_1A_rpc_adjusted = test[i]$D1_1A_rpc_adjusted)
fast_peaks_dt <- cbind(fast_peaks_dt, D1_2A_rpc_adjusted = test[i]$D1_2A_rpc_adjusted)
fast_peaks_dt <- cbind(fast_peaks_dt, D6_2A_rpc_adjusted = test[i]$D6_2A_rpc_adjusted)
fast_peaks_dt <- cbind(fast_peaks_dt, D12_1A_rpc_adjusted = test[i]$D12_1A_rpc_adjusted)
fast_peaks_dt <- cbind(fast_peaks_dt, D12_2A_rpc_adjusted = test[i]$D12_2A_rpc_adjusted)
fast_peaks_dt[, D1_1A_norm := D1_1A / D1_1A_rpc_adjusted]
fast_peaks_dt[, D1_2A_norm := D1_2A / D1_2A_rpc_adjusted]
fast_peaks_dt[, D6_2A_norm := N2_D6_2A / D6_2A_rpc_adjusted]
fast_peaks_dt[, D12_1A_norm := D12_1A / D12_1A_rpc_adjusted]
fast_peaks_dt[, D12_2A_norm := D12_2A / D12_2A_rpc_adjusted]
fast_peaks_dt[, D1_norm := (D1_1A_norm + D1_2A_norm) / 2]
fast_peaks_dt[, D12_norm := (D12_1A_norm + D12_2A_norm) / 2]
fast_peaks_dt[, position_norm := position / length]

plot <- ggplot(fast_peaks_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                               D6_2A_rpc_adjusted >= 1 & 
                               D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 &
                               peak > 25 & (peak - length) < -25,
                             .(adjusted, D1 = movingAverage(D1_norm, n=2, center=T),
                               D6 = movingAverage(D6_2A_norm, n=2, center=T),
                               D12 = movingAverage(D12_norm, n=2, center=T))]) +
  stat_summary(aes(adjusted, D6), fun.y = "mean", geom = "line", size=1.25, color = '#F28CC1') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=1.25, color = '#E7298A') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_color_manual(labels = c("Day 1","Day 6", "Day 12"), values = c("gray40", "#F28CC1", "#E7298A"), name = "") +
  scale_x_continuous(expand = expansion(), limits = c(-25, 25))
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/MetageneFastPeaks_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

library("Logolas")
logomaker(fast_peaks[D1_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)

bg_noW <- bg[c(1:18,20)]
logomaker(motif_dt[count >= 100 & ratio_12_1 < 1 & D1_pause > 2]$motif, type = "EDLogo", bg = bg_noW, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)

logomaker(motif_dt[count >= 100 & ratio_6_1 < 1 & D1_pause > 2]$motif, type = "EDLogo", bg = bg_noW, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_6_1 > 1 & D6_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)

logomaker(motif_dt[count >= 100 & ratio_12_6 > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg_noW, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_12_6 < 1 & D6_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)

weblogo(fast_peaks[D1_pause > 10]$motif3, open = TRUE, file.out = "fast_probability.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
weblogo(fast_peaks[D1_pause > 10]$motif3, open = TRUE, file.out = "fast_IC.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        composition = c('A'=6.325907,'C'=2.076066,'D'=5.304005,'E'=6.526124,'F'=4.787855,'G'=5.344524,'H'=2.281571,'I'=6.213785,'K'=6.339429,'L'=8.643033,'M'=2.648252,'N'=4.877476,'P'=4.892958,'Q'=4.081622,'R'=5.137781,'S'=8.080618,'T'=5.891986,'V'=6.229459,'W'=1.107328,'Y'=3.210222),
        color.scheme = 'chemistry', errorbars = FALSE, yaxis = 0.05)



### Proline PP and PPP ###
PPP_dt <- readRDS("doc/Archive_Misc/PPP_dt.rds")
PP_dt <- readRDS("doc/Archive_Misc/PP_dt.rds")
plot <- ggplot(data = PP_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                              D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1,
                            .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=3, center=T),
                              D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7298A') +
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size = 1.25, color = '#E7298A') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size = 1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from PP (codons)") +
  theme(legend.position = "none", legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Worm_PP.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


