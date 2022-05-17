library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')

### Figure 1

# Fig 1b: Global pausing
N2_dtA <- readRDS("doc/N2_dtA.rds")
plot <- ggplot(N2_dtA[D1_1A_rpc >= 1 & D12_1A_rpc >= 1 & 
                D1_2A_rpc >= 1 & D12_2A_rpc >= 1 & 
                N2_D6_2A_rpc >= 1 &
                D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
                N2_D6_2A_sum >= 64 &  
                D1_1A_coverage > 0.7 & D1_2A_coverage > 0.7 &
                D12_1A_coverage > 0.7 & D12_2A_coverage > 0.7 & 
                N2_D6_2A_coverage > 0.7 &
                position > 20 & stopdist < -20], 
       aes(D1_pause)) + scale_x_log10() + coord_cartesian(xlim = c(0.05,20)) +
  stat_ecdf(aes(N2_D6_2A_pause), color = "#F28CC1", geom = "step", size = 1.25) +
  stat_ecdf(aes(D12_pause), color = "#E7298A", geom = "step", size = 1.25) +
  stat_ecdf(geom = "step", size = 1.25, color = "gray40") +
  scale_color_manual(labels = c("Day 1","Day 6", "Day 12"), values = c("gray40", "#F28CC1", "#E7298A"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Cumulative fraction", x = "Pause score") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/globalPausing_worm.tiff", plot, width = 6, height = 4, dpi = 300)

# Fig 1c: Volcano of age-dependent stalls
N2_fishers <- readRDS("doc/N2_fishers.rds")
stalling_peaks_all <- readRDS("stalling_peaks_all.rds")
plot <- ggplot(N2_fishers[odds < Inf & odds > 0 & padj != 0]) + 
  geom_point(aes(log2(odds), -log10(padj)), color = "gray80", alpha = 0.5, size = 2) +
  #geom_point(data = stalling_peaks_all[D12_pause < N2_D6_2A_pause], aes(log2(odds), -log10(padj), color = "#9E9AC8"), alpha = 0.5, size = 2) +
  geom_point(data = stalling_peaks_all[D12_pause > N2_D6_2A_pause & padj != 0], aes(log2(odds), -log10(padj), color = "#E7298A"), alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous() + 
  scale_color_manual(labels = c(""), values = c("#E7298A"), name = "")
  #scale_color_manual(labels = c("Pausing","Age-dep. pausing"), values = c("#9E9AC8", "#FF7F00"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D12/D1)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_worms.tiff", plot, width = 6, height = 4, dpi = 300)

# Fig 1d: Metagene of age-dependent pausing sites
stalling_peaks_ageDep_dt <- readRDS("stalling_peaks_ageDep_dt.rds")
plot <- ggplot(stalling_peaks_ageDep_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
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
ggsave("/Users/KevinStein/Desktop/stallingMetagene_worms_ageDep.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure S1

# FigS1a: Polysome profiles
polysome <- read.csv("/Users/KevinStein/Desktop/Manuscript/Figures/Plots/FigS1/polysomeProfilesWorms.csv", stringsAsFactors = T, header = T)

plot <- ggplot(polysome, aes(Position, D1_1)) + geom_line(size = 1.25, color = "gray30") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/D1_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(polysome, aes(Position, D12_1)) + geom_line(size = 1.25, color = "#E31A1C") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/D12_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(polysome, aes(Position, (N2_D6_1+1))) + geom_line(size = 1.25, color = "#1F78B4") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/D6_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

### Figure S2

N2_expression_deseq <- readRDS("doc/N2_expression_deseq.rds")

# FigS2a: Replicates
plot <- ggplot(N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64], aes(D1_1A_tpm, D1_2A_tpm)) +
  geom_point(aes(D1_1A_tpm, D1_2A_tpm), fill = "gray30", color = "gray30", alpha = 0.5, size = 2) + 
  scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
  scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
  annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Day 1, rep 2 (log10 TPM)", x = "Day 1, rep 1 (log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
cor(N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64]$D1_1A_tpm,
    N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64]$D1_2A_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/D1_correlation_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# ggplot(N2_expression[D1_1A_sum > 64 & N2_D1_2A_sum > 64], aes(D1_1A_tpm, N2_D1_2A_tpm)) + geom_point() +
#   scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5")) + 
#   scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5")) + 
#   annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))

plot <- ggplot(N2_expression_deseq[D12_1A_sum > 64 & D12_2A_sum > 64], aes(D12_1A_tpm, D12_2A_tpm)) +
  geom_point(aes(D12_1A_tpm, D12_2A_tpm), fill = "gray30", color = "gray30", alpha = 0.5, size = 2) + 
  scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
  scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
  annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Day 12, rep 2 (log10 TPM)", x = "Day 12, rep 1 (log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
cor(N2_expression_deseq[D12_1A_sum > 64 & D12_2A_sum > 64]$D12_1A_tpm,
    N2_expression_deseq[D12_1A_sum > 64 & D12_2A_sum > 64]$D12_2A_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/D12_correlation_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# FigS2c: Gene-level expression
plot <- ggplot(N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 0]) + 
  geom_point(aes(log2FoldChange, -log10(padj)), fill = "gray75", color = "gray75", alpha = 0.5, size = 2) + 
  geom_point(data = N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 1], aes(log2FoldChange, -log10(padj)), color = factor(group), fill = "#E31A1C", alpha = 0.5, size = 2) + 
  geom_point(data = N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 2], aes(log2FoldChange, -log10(padj), color = factor(group)), fill = "#377EB8", alpha = 0.5, size = 2) + 
  scale_color_manual(limits = c("1","2"), labels = c(">2 fold up",">2 fold down"), values = c("#9970AB", "#5AAE61"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "log2 (fold change)") +
  theme(legend.position = c(.01,1.1), legend.justification = c("left", "top"), legend.title=element_text(size=16), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Deseq.gene.expression_volcano_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# plot <- ggplot(N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 0]) + 
#   geom_point(aes(D1_tpm, D12_tpm), fill = "gray75", color = "gray75", alpha = 0.5, size = 2) + 
#   geom_point(data = N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 1], aes(D1_tpm, D12_tpm, color = factor(group)), fill = "#E31A1C", alpha = 0.5, size = 2) + 
#   geom_point(data = N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 2], aes(D1_tpm, D12_tpm, color = factor(group)), fill = "#377EB8", alpha = 0.5, size = 2) + 
#   scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
#   scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
#   annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm")) +
#   scale_color_manual(limits = c("1","2"), labels = c(">2 fold up",">2 fold down"), values = c("#9970AB", "#5AAE61"), name = "Pearson's r = 0.83")
# plot <- plot + theme_classic(20) + labs(y = "Day 12 (log10 TPM)", x = "Day 1 (log10 TPM)") +
#   theme(legend.position = c(.53,.4), legend.justification = c("left", "top"), legend.title=element_text(size=16), legend.background = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
# cor(N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64]$D1_tpm,
#     N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64]$D12_tpm, method = "pearson")
# ggsave("/Users/KevinStein/Desktop/Deseq.gene.expression_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# FigS2e: Gene Ontology
GO <- read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/GO_translatome.csv", stringsAsFactors = T, header = T)
GO <- as.data.table(GO)
plot <- ggplot(GO[Organism == "worm"], aes(x = reorder(Description, -foldchange), y = foldchange, fill = factor(Group))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("down","up"), values = c("#5AAE61", "#9970AB"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_BP_translatomeWorm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure 4 ###

# Fig 4a: Sequence logo
stalling_peaks_ageDep <- readRDS("stalling_peaks_ageDep.rds")
logomaker(stalling_peaks_ageDep[D12_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(stalling_peaks_all_ageDep[D12_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
#logomaker(stalling_peaks_ageDep[D12_pause > 10]$motif3, type = "Logo", bg = bg, color_seed = 6)
# temp <- stalling_peaks_ageDep_dt[position == peak_start]
# logomaker(temp[D12_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)

# Fig 4b: Pausing at xRR
RR_dt <- readRDS("RR_dt.rds")
plot <- ggplot(data = RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
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
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = "none", legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Worm_RR_smooth3.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                    D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 & 
                    adjusted == 0]$D1_2A_norm,
            RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                    D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 & 
                    adjusted == 0]$D12_2A_norm, alternative = 't')

# Fig 4c: Pausing at polybasic KR6
KR6of6_dt <- readRDS("doc/KR6of6_dt.rds")

plot <- ggplot(data = KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                                  D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1,
                                .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=3, center=T),
                                  D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7298A') +
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size = 1.25, color = '#E7298A') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size = 1.25, color = 'gray40') +
  coord_cartesian(xlim = c(-25, 25), ylim = c(0.8,1.35), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = "none", legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Worm_KR6of6.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Fig 4d: tag-342 example
plot <- ggplot(data=N2_dtA[orf == "B0464.8", .(position, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=5, center=T), 
                                               D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=5, center=T))]) + 
  geom_rect(xmin = 243, xmax = 283, ymin = -Inf, ymax = Inf, fill = "#F28CC1", alpha = 0.1) +
  geom_vline(xintercept = 263, color = "gray50", linetype = "longdash", size = 1) +
  geom_line(aes(position, D12), color = "#E7298A", size = 1.25) +
  geom_line(aes(position, D1), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expansion())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/B0464.8_KR5of5_aa263.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=N2_dtA[orf == "B0464.8", .(position, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=2, center=T), 
                                               D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=2, center=T))]) + 
  geom_vline(xintercept = 263, color = "gray50", linetype = "longdash", size = 1) +
  geom_line(aes(position, D12), color = "#E31A1C", size = 1.25) +
  geom_line(aes(position, D1), color = "gray30", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(243,283)) + coord_cartesian(ylim = c(0,5.5))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/B0464.8_KR5of5_aa263_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S5g: K12C11.6 example, W in P-site
plot <- ggplot(data=N2_dtA[orf == "K12C11.6", .(position, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=3, center=T), 
                                               D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=3, center=T))]) + 
  #geom_vline(xintercept = 62, color = "gray50", linetype = "longdash", size = 1) +
  geom_line(aes(position, D12), color = "#E7298A", size = 1.25) +
  geom_line(aes(position, D1), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expansion())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/K12C11.6_Waa62.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig 4e: Association of aggregation and pausing
temp <- data.table(Category = c("No", "Stall"), 
                   FractionAggregation = c(((length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                                               length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)) /
                                              (length(hartl_proteome_orfs$WBgene) - 4)),
                                           (length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)) /
                                             length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene)))
plot <- ggplot(temp, aes(x = Category, y=FractionAggregation, fill = Category)) + geom_col(color = "black", size = 0.5) + scale_y_continuous(limits = c(0,0.25), expand = c(0.0008,0.0008)) +
  scale_fill_manual(limits = c("No","Stall"), values = c("#80B1D3", "#FDB462"), name = "")
plot <- plot + theme_classic(18) + labs(y = "Fraction of aggregated\nproteins in dataset", x = "") +
  theme(legend.position = "none", legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Stalling_HartlAggregation_p5.6e-22.pdf", plot, width = 3.2, height = 4, dpi = 300, useDingbats = F)

# Fig 4f: GO of genes with pausing and aggregation
plot <- ggplot(GO_Hartl_subset, aes(x = reorder(Description, div), y = Fold.Enrichment, fill = factor(Group))) + 
  geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("all","stall"), values = c("#80B1D3", "#FDB462"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = c(0.85,0.15), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_Hartl_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Fig 4g: wars-1
plot <- ggplot(data=N2_dtA[gene == "wars-1", .(position, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=5, center=T), 
                                               D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=5, center=T))]) + 
  geom_rect(xmin = 146, xmax = 186, ymin = -Inf, ymax = Inf, fill = "#F28CC1", alpha = 0.1) +
  geom_line(aes(position, D12), color = "#E7298A", size = 1.25) +
  geom_line(aes(position, D1), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/wars-1_aa166.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=N2_dtA[gene == "wars-1", .(position, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=5, center=T), 
                                               D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=5, center=T))]) + 
  geom_vline(xintercept = 166, color = "gray50", linetype = "longdash", size = 1) +
  geom_line(aes(position, D12), color = "#E31A1C", size = 1.25) +
  geom_line(aes(position, D1), color = "gray30", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(146,186))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/wars-1_aa166_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure S5 ###

# FigS5a-c: Sequence logos
motif_dt <- readRDS("motif_dt.rds")
logomaker(motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)
#logomaker(motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2]$motif, type = "Logo", bg = bg, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_6_1 > 1 & D6_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)
#logomaker(motif_dt[count >= 100 & ratio_6_1 > 1 & D6_pause > 2]$motif, type = "Logo", bg = bg, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_12_6 > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg_noW, color_seed = 6)
#logomaker(motif_dt[count >= 100 & ratio_12_6 > 1 & D12_pause > 2]$motif, type = "Logo", bg = bg_noW, color_seed = 6)

# FigS5d: Residue frequency
plot <- ggplot(aa_freq[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-0.9,0.6))
plot <- plot + theme_classic(20) + labs(y = "A-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/A-siteFreq_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aa_freqP[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-0.9,0.6))
plot <- plot + theme_classic(20) + labs(y = "P-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/P-siteFreq_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aa_freqE[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-0.9,0.6))
plot <- plot + theme_classic(20) + labs(y = "E-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/E-siteFreq_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

test <- matrix(c(aa_freq[residue == "V"]$stalls,
                 aa_freq[residue == "V"]$proteome,
                 (sum(aa_freq$stalls) - aa_freq[residue == "V"]$stalls),
                 (sum(aa_freq$proteome) - aa_freq[residue == "V"]$proteome)), nrow = 2)
fisher.test(test) # R: p-value = 0.004836
test <- matrix(c(aa_freqP[residue == "F"]$stalls,
                 aa_freqP[residue == "F"]$proteome,
                 (sum(aa_freqP$stalls) - aa_freqP[residue == "F"]$stalls),
                 (sum(aa_freqP$proteome) - aa_freqP[residue == "F"]$proteome)), nrow = 2)
fisher.test(test) # R: p-value = 0.01927
test <- matrix(c(aa_freqE[residue == "E"]$stalls,
                 aa_freqE[residue == "E"]$proteome,
                 (sum(aa_freqE$stalls) - aa_freqE[residue == "E"]$stalls),
                 (sum(aa_freqE$proteome) - aa_freqE[residue == "E"]$proteome)), nrow = 2)
fisher.test(test) # R: p-value = 0.04878

# FigS5e: Codon frequency
plot <- ggplot(codon_freq[codon != 'TAA' & codon != 'TAG' & codon != 'TGA'], aes(x = codon, y = log2(div), fill = residue)) + geom_col() +
  geom_hline(yintercept = 0, size = 0.2, color = 'black')
plot <- plot + theme_minimal(20, base_line_size = 0.5) + labs(y = "A-site frequency, log2", x = "Codon") +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 5, color = 'black'), axis.ticks.x = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 16, color = "black"), axis.ticks.y = element_line(color = "black"))
ggsave("/Users/KevinStein/Desktop/codonfreq_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# FigS5f: Polybasic analysis
KR4of4_dt <- readRDS("doc/KR4of4_dt.rds")
KR5of5_dt <- readRDS("doc/KR5of5_dt.rds")

plot <- ggplot(data = KR4of4_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                                  D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1,
                                .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=3, center=T),
                                  D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7298A') +
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size = 1.25, color = '#E7298A') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size = 1.25, color = 'gray40') +
  #coord_cartesian(xlim = c(-25, 25), ylim = c(0.8,1.35), expand = expand_scale()) +
  coord_cartesian(xlim = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = "none", legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Worm_KR4of4.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                                  D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1,
                                .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=3, center=T),
                                  D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7298A') +
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size = 1.25, color = '#E7298A') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size = 1.25, color = 'gray40') +
  #coord_cartesian(xlim = c(-25, 25), ylim = c(0.8,1.35), expand = expand_scale()) +
  coord_cartesian(xlim = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = "none", legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Worm_KR5of5.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Stats
# KR4of4_dt[, D1_norm := (D1_1A_norm + D1_2A_norm) / 2]
# KR4of4_dt[, D12_norm := (D12_1A_norm + D12_2A_norm) / 2]
# KR5of5_dt[, D1_norm := (D1_1A_norm + D1_2A_norm) / 2]
# KR5of5_dt[, D12_norm := (D12_1A_norm + D12_2A_norm) / 2]
# KR6of6_dt[, D1_norm := (D1_1A_norm + D1_2A_norm) / 2]
# KR6of6_dt[, D12_norm := (D12_1A_norm + D12_2A_norm) / 2]
# RR_dt[, D1_norm := (D1_1A_norm + D1_2A_norm) / 2]
# RR_dt[, D12_norm := (D12_1A_norm + D12_2A_norm) / 2]
# 
# plot <- ggplot(data = KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                                   D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 3 & adjusted < 6], 
#                aes("6_1", D1_norm)) + 
#   geom_boxplot(fill = "gray50", notch = T) +
#   geom_boxplot(data = KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                                   D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 3 & adjusted < 6], 
#                aes("6_2", D12_norm), fill = "#E7298A", notch = T) +
#   coord_capped_cart(left = "both", ylim = c(0,4)) +
#   scale_x_discrete(limits = c("6_1", "6_2"), labels = c("Day 1", "Day 12"))
# plot <- plot + theme_classic(18) + labs(y = "Norm. ribosome occupancy", x = "") +
#   theme(legend.position = c(0.9,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
#         axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
#         plot.background = element_blank())
# ggsave("/Users/KevinStein/Desktop/KR6of6_p0.002601.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)
# 
# ggplot(data = KR4of4_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                           D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 1 & adjusted < 4], 
#        aes("4_1", D1_norm)) + 
#   geom_boxplot() +
#   geom_boxplot(data = KR4of4_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                                   D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 1 & adjusted < 4], 
#                aes("4_2", D12_norm)) +
#   geom_boxplot(data = KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                                   D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 2 & adjusted < 5], 
#                aes("5_1", D1_norm)) +
#   geom_boxplot(data = KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                                   D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 2 & adjusted < 5], 
#                aes("5_2", D12_norm)) +
#   geom_boxplot(data = KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                                   D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 3 & adjusted < 6], 
#                aes("6_1", D1_norm)) +
#   geom_boxplot(data = KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                                   D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 3 & adjusted < 6], 
#                aes("6_2", D12_norm)) +
#   coord_cartesian(ylim = c(0,4))
# 
# wilcox.test(KR4of4_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 1 & adjusted < 4]$D1_norm,
#             KR4of4_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 1 & adjusted < 4]$D12_norm, alternative = "t")
# wilcox.test(KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 2 & adjusted < 5]$D1_norm,
#             KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 2 & adjusted < 5]$D12_norm, alternative = "t")
# wilcox.test(KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 3 & adjusted < 6]$D1_norm,
#             KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted > 3 & adjusted < 6]$D12_norm, alternative = "t")
# wilcox.test(KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted == 6]$D1_norm,
#             KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted == 6]$D12_norm, alternative = "t")
# wilcox.test(RR_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                     D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted == 0]$D1_norm,
#             KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
#                         D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & adjusted == 0]$D12_norm, alternative = "t")

# FigS5g: GO of genes with pausing
GO <- read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/GO_stalling.csv", stringsAsFactors = T, header = T)
GO <- as.data.table(GO)
plot <- ggplot(GO[Organism == "worms"], aes(x = reorder(Description, Fold.Enrichment), y = Fold.Enrichment, fill = factor(Category))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("BP","CC"), values = c("#A6CEE3", "#B2DF8A"), name = "")
plot <- plot + theme_classic(12) + labs(y = "Fold enrichment", x = "") +
  theme(legend.position = c(0.8,0.2), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_stallingWorm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure S6 ###

# FigS6a: Association of pausing and aggregation
temp <- data.table(Category = c("No", "Stall"), 
                   FractionAggregation = c(((length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                                               length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)) /
                                              (length(hartl_proteome_orfs$WBgene) - 4)),
                                           (length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)) /
                                             length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene)))
plot <- ggplot(temp, aes(x = Category, y=FractionAggregation, fill = Category)) + geom_col(color = "black", size = 0.5) + scale_y_continuous(limits = c(0,0.15), expand = c(0.0005,0.0005)) +
  scale_fill_manual(limits = c("No","Stall"), values = c("#80B1D3", "#FDB462"), name = "")
plot <- plot + theme_classic(18) + labs(y = "Fraction of aggregated\nproteins in dataset", x = "") +
  theme(legend.position = "none", legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Stalling_kenyonAggregation_p5.6e-22.pdf", plot, width = 3.2, height = 4, dpi = 300, useDingbats = F)

# FigS6b: Enrichment of basic residues
plot <- ggplot(agg_polybasic_dt[datasource == "Hartl"], aes(x = residues, y = odds, fill = stall)) + 
  geom_col(color = "black", size = 0.5, position = "dodge") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(limits = c("a","y"), labels = c("all", "yes"), values = c("gray50", "#00CED1"), name = "Pause site") +
  scale_y_continuous(expand = c(0.0008,0.0008)) +
  scale_x_discrete(labels = c("4 Lys/Arg", "2 Arg"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio (relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.1,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_HartlAggregation.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)

plot <- ggplot(agg_polybasic_dt[datasource == "Kenyon"], aes(x = residues, y = odds, fill = stall)) + 
  geom_col(color = "black", size = 0.5, position = "dodge") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(limits = c("a","y"), labels = c("all", "yes"), values = c("gray50", "#00CED1"), name = "Pause site") +
  scale_y_continuous(expand = c(0.0008,0.0008)) +
  scale_x_discrete(labels = c("4 Lys/Arg", "2 Arg"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio (relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.1,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_KenyonAggregation.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)

# FigS5c: GO of genes with pausing and aggregation
plot <- ggplot(GO_Kenyon_subset, aes(x = reorder(Description, div), y = Fold.Enrichment, fill = factor(Group))) + 
  geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("all","stall"), values = c("#FB9A99", "#E31A1C"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = c(0.85,0.15), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_Kenyon_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# FigS5d: Aggregation of tRNA synthetases
plot <- ggplot(tRNA_aggregation1, aes(age,log2.ratio.lh, color = gene, group = gene)) + geom_line(size = 1) + geom_point(size = 2)
plot <- plot + theme_classic(18) + labs(y = "Aggregate abundance, log2", x = "Age (days)") +
  theme(legend.position = c(0.75,0.4), legend.text = element_text(size = 12), legend.title = element_blank(), legend.background = element_blank(),
        axis.text = element_text(color = "black", size = 14),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/tRNA_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Composition of Zpeaks
# plot <- ggplot(freq[residue != "X"], aes(residue,log2(freq), fill = age)) + 
#   geom_col(position = "dodge") + 
#   scale_fill_manual(limits = c("D1", "D6", "D12"), 
#                     labels = c("Day 1", "Day 6", "Day 12"),
#                     values = c("#1F78B4", "#999999", "#E31A1C"), name = "")
# plot <- plot + theme_classic(20) + labs(y = "log2(frequency)", x = "Residue") +
#   theme(legend.position = "none",
#         axis.text = element_text(size = 16, color = "black"))
# ggsave("/Users/KevinStein/Desktop/frequency_Zpeaks_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Metagene from start
# plot <- ggplot(data = N2_dtA[D1_1A_rpc >= 1 & D12_1A_rpc >= 1 & 
#                        D1_2A_rpc >= 1 & D12_2A_rpc >= 1 & 
#                        N2_D6_2A_rpc >= 1 &
#                        D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
#                        D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
#                        N2_D6_2A_sum >= 64 & length >= 300,
#                      .(position, D1 = movingAverage(D1_pause, n=2, center=T),
#                        D6 = movingAverage(N2_D6_2A_pause, n=2, center=T),
#                        D12 = movingAverage(D12_pause, n=2, center=T))]) + xlim(-7, 200) + coord_cartesian(ylim = c(0.9,1.4), expand = expand_scale()) +
#   stat_summary(aes(position, D6, color = "#999999"), fun.y = "mean", geom = "line", size = 1.25) +
#   stat_summary(aes(position, D1, color = "#1F78B4"), fun.y = "mean", geom = "line", size = 1.25) +
#   stat_summary(aes(position, D12, color = "#E31A1C"), fun.y = "mean", geom = "line", size = 1.25) +
#   scale_color_manual(labels = c("Day 1","Day 6", "Day 12"), values = c("#1F78B4", "#999999", "#E31A1C"), name = "")
# plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Codon position") +
#   theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
#         panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
# ggsave("/Users/KevinStein/Desktop/startMetagene_worm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Metagene from stop
# ggplot(data = N2_dtA[D1_1A_rpc >= 1 & D12_1A_rpc >= 1 & 
#                                D1_2A_rpc >= 1 & D12_2A_rpc >= 1 & 
#                                N2_D6_2A_rpc >= 1 &
#                                D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
#                                D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
#                                N2_D6_2A_sum >= 64 & length >= 300,
#                              .(stopdist, D1 = movingAverage(D1_pause, n=2, center=T),
#                                D6 = movingAverage(N2_D6_2A_pause, n=2, center=T),
#                                D12 = movingAverage(D12_pause, n=2, center=T))]) + xlim(-200, 7) + 
#   stat_summary(aes(stopdist, D6, color = "#999999"), fun.y = "mean", geom = "line", size = 1.25) +
#   stat_summary(aes(stopdist, D1, color = "#1F78B4"), fun.y = "mean", geom = "line", size = 1.25) +
#   stat_summary(aes(stopdist, D12, color = "#E31A1C"), fun.y = "mean", geom = "line", size = 1.25) +
#   scale_color_manual(labels = c("Day 1","Day 6", "Day 12"), values = c("#1F78B4", "#999999", "#E31A1C"), name = "")
# plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Codon position") +
#   theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
#         panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
# ggsave("/Users/KevinStein/Desktop/stopMetagene_worm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# temp <- D12_asymmetry_max[i.position > 125 & (i.position - length) < -35 & D12_asymmetry < 1 & D1_asymmetry > 1, c(1)]
# for(i in temp$orf) {
#   filename <- paste("/Users/KevinStein/Desktop/test/", i, ".pdf", sep = "")
#   G <- ggplot(data=N2_dtA[orf == i, .(position, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=5, center=T),
#                                             D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=5, center=T))]) +
#     geom_vline(data = D12_asymmetry_max[orf == i], xintercept = D12_asymmetry_max[orf == i]$position, color = "blue", linetype = "dashed") +
#     geom_line(aes(position, D12), color = "red") +
#     geom_line(aes(position, D1), color = "black")
#   ggsave(filename, G, width = 6, height = 4, dpi = 300)
# }


