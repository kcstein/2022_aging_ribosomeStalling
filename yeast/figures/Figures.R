### Figures ###
library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
sc_dtA <- readRDS("doc/sc_dtA.rds")
sc_fishers <- readRDS("doc/sc_fishers.rds")

# Fig 1b: Global stalling
plot <- ggplot(sc_dtA[WT_D0_1A_rpc >= 1 & WT_D0_2A_rpc >= 1 & 
                        WT_D2_1A_rpc >= 1 & WT_D2_2A_rpc >= 1 &
                        WT_D4_1A_rpc >= 1 & WT_D4_2A_rpc >= 1 &
                        WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                        WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
                        WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 &
                        position > 20 & stopdist < -20], 
               aes(WT_D0_pause)) + scale_x_log10() + coord_cartesian(xlim = c(0.05,20)) +
  stat_ecdf(aes(WT_D2_pause), color = "#F2CF8C", geom = "step", size = 1.25) +
  stat_ecdf(aes(WT_D4_pause), color = "#E7A427", geom = "step", size = 1.25) +
  stat_ecdf(geom = "step", size = 1.25, color = "gray40") +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray40", "#F2CF8C", "#E7A427"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Cumulative fraction", x = "Pause score") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/globalPausing_yeast.tiff", plot, width = 6, height = 4, dpi = 300)

# Fig 1c: Volcano of age-dependent stalls
sc_fishers <- readRDS("doc/sc_fishers.rds")
WT_peaks_all <- readRDS("doc/WT_peaks_all.rds")
sch9_peaks_all <- readRDS("sch9_peaks_all.rds")
plot <- ggplot(sc_fishers[WT_odds < Inf & WT_odds > 0 & WT_padj != 0]) + 
  geom_point(aes(log2(WT_odds), -log10(WT_padj)), color = "gray80", alpha = 0.5, size = 2) +
  #geom_point(data = WT_peaks_all[WT_D4_pause < WT_D2_pause], aes(log2(WT_odds), -log10(WT_padj), color = "#9E9AC8"), alpha = 0.5, size = 2) +
  geom_point(data = WT_peaks_all[WT_D4_pause > WT_D2_pause & WT_padj != 0], aes(log2(WT_odds), -log10(WT_padj), color = "#E7A427"), alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous(limits = c(0,200)) + 
  scale_color_manual(labels = c(""), values = c("#E7A427"), name = "")
  #scale_color_manual(labels = c("Pausing","Age-dep. pausing"), values = c("#9E9AC8", "#FF7F00"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D4/D0)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_yeast.tiff", plot, width = 6, height = 4, dpi = 300)

# plot <- ggplot(sc_fishers[sch9_odds < Inf & sch9_odds > 0]) + 
#   geom_point(aes(log2(sch9_odds), -log10(sch9_padj)), color = "gray80", alpha = 0.5, size = 2) +
#   #geom_point(data = WT_peaks_all[WT_D4_pause < WT_D2_pause], aes(log2(WT_odds), -log10(WT_padj), color = "#9E9AC8"), alpha = 0.5, size = 2) +
#   geom_point(data = sch9_peaks_all[sch9_D4_pause > sch9_D2_pause], aes(log2(sch9_odds), -log10(sch9_padj), color = "#E31A1C"), alpha = 0.6, size = 2) +
#   scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(0,200)) + 
#   scale_color_manual(labels = c("Age-dep. pausing"), values = c("#E31A1C"), name = "")
# #scale_color_manual(labels = c("Pausing","Age-dep. pausing"), values = c("#9E9AC8", "#FF7F00"), name = "")
# plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D4/D0)") +
#   theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
#         panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
# ggsave("/Users/KevinStein/Desktop/Volcano_stalling_sch9.tiff", plot, width = 6, height = 4, dpi = 300)


# Fig 1d: Metagene of age-dependent pausing sites
WT_stalling_peaks_ageDep_dt <- readRDS("WT_stalling_peaks_ageDep_dt.rds")
plot <- ggplot(data = WT_stalling_peaks_ageDep_dt[WT_D0_1A_rpc_adjusted >= 1 & WT_D0_2A_rpc_adjusted >= 1 & 
                                            WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                                            WT_D4_1A_rpc_adjusted >= 1 & WT_D4_2A_rpc_adjusted >= 1 &
                                            peak > 25 & (peak - length) < -25,
                                          .(adjusted, WT0 = movingAverage(WT_D0_norm, n=1, center=T),
                                            WT2 = movingAverage(WT_D2_norm, n=1, center=T),
                                            WT4 = movingAverage(WT_D4_norm, n=1, center=T))]) +
  # stat_summary(aes(adjusted, WT2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
  #              fun.args=list(conf.int=0.5), fill = '#F2CF8C') +
  # stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
  #              fun.args=list(conf.int=0.5), fill = '#E7A427') +
  # stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
  #              fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = '#F2CF8C') +  
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray40", "#F2CF8C", "#E7A427"), name = "") +
  scale_x_continuous(expand = expand_scale(), limits = c(-25,25))
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_yeast_ageDep.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure S1

# FigS1b: Polysome profiles
polysome <- read.csv("/Users/KevinStein/Desktop/Manuscript/1st Submission/Data/Plots/FigS1/polysomeProfilesYeast.csv", stringsAsFactors = T, header = T)

plot <- ggplot(polysome, aes(Position, (WT_D0_2 + 1))) + geom_line(size = 1.25, color = "gray30") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/yeastD0_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(polysome, aes(Position, (WT_D2_2 + 1))) + geom_line(size = 1.25, color = "#1F78B4") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/yeastD2_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(polysome, aes(Position, (WT_D4_2 + 1))) + geom_line(size = 1.25, color = "#E31A1C") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/yeastD4_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure S2

WT_expression <- readRDS("doc/WT_expression.rds")

# Fig S2b: Replicates
WT_expression <- readRDS("doc/WT_expression.rds")
plot <- ggplot(WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64], aes(WT_D0_1A_tpm, WT_D0_2A_tpm)) +
  geom_point(aes(WT_D0_1A_tpm, WT_D0_2A_tpm), fill = "gray30", color = "gray30", alpha = 0.5, size = 2) + 
  scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
  scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
  annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Day 0, rep 2 (log10 TPM)", x = "Day 0, rep 1 (log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
cor(WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64]$WT_D0_1A_tpm,
    WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64]$WT_D0_2A_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/D0_correlation_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(WT_expression[WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64], aes(WT_D4_1A_tpm, WT_D4_2A_tpm)) +
  geom_point(aes(WT_D4_1A_tpm, WT_D4_2A_tpm), fill = "gray30", color = "gray30", alpha = 0.5, size = 2) + 
  scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
  scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
  annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Day 4, rep 2 (log10 TPM)", x = "Day 4, rep 1 (log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
cor(WT_expression[WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64]$WT_D4_1A_tpm,
    WT_expression[WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64]$WT_D4_2A_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/D4_correlation0.96_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# ggplot(expression[sch9_D0_1A_sum > 64 & sch9_D0_2A_sum > 64], aes(sch9_D0_1A_tpm, sch9_D0_2A_tpm)) + geom_point() +
#   scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
# cor(expression[sch9_D0_1A_sum > 64 & sch9_D0_2A_sum > 64]$sch9_D0_1A_tpm,
#     expression[sch9_D0_1A_sum > 64 & sch9_D0_2A_sum > 64]$sch9_D0_2A_tpm, method = "pearson")
# 
# ggplot(expression[sch9_D4_1A_sum > 64 & sch9_D4_2A_sum > 64], aes(sch9_D4_1A_tpm, sch9_D4_2A_tpm)) + geom_point() +
#   scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
# cor(expression[sch9_D4_1A_sum > 64 & sch9_D4_2A_sum > 64]$sch9_D4_1A_tpm,
#     expression[sch9_D4_1A_sum > 64 & sch9_D4_2A_sum > 64]$sch9_D4_2A_tpm, method = "pearson")

# FigS1d: Gene-level expression
plot <- ggplot(WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 0]) + 
  geom_point(aes(log2FoldChange, -log10(padj)), fill = "gray75", color = "gray75", alpha = 0.5, size = 2) + 
  geom_point(data = WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 1], aes(log2FoldChange, -log10(padj), color = factor(group)), fill = "#E31A1C", alpha = 0.5, size = 2) + 
  geom_point(data = WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 2], aes(log2FoldChange, -log10(padj), color = factor(group)), fill = "#377EB8", alpha = 0.5, size = 2) + 
  scale_color_manual(limits = c("1","2"), labels = c(">2 fold up",">2 fold down"), values = c("#9970AB", "#5AAE61"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "log2 (fold change)") +
  theme(legend.position = c(.01,1.1), legend.justification = c("left", "top"), legend.title=element_text(size=16), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Deseq.gene.expression_volcano_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# plot <- ggplot(WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 0]) + 
#   geom_point(aes(D0_tpm, D4_tpm), fill = "gray75", color = "gray75", alpha = 0.5, size = 2) + 
#   geom_point(data = WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 1], aes(D0_tpm, D4_tpm, color = factor(group)), fill = "#E31A1C", alpha = 0.5, size = 2) + 
#   geom_point(data = WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 2], aes(D0_tpm, D4_tpm, color = factor(group)), fill = "#377EB8", alpha = 0.5, size = 2) + 
#   scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
#   scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0", "1","2","3", "4","5"), expand = expand_scale()) + 
#   annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm")) +
#   scale_color_manual(limits = c("1","2"), labels = c(">2 fold up",">2 fold down"), values = c("#E31A1C", "#1F78B4"), name = "Pearson's r = 0.39")
# plot <- plot + theme_classic(20) + labs(y = "Day 4 (log10 TPM)", x = "Day 0 (log10 TPM)") +
#   theme(legend.position = c(.53,.4), legend.justification = c("left", "top"), legend.title=element_text(size=16), legend.background = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
# cor(WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64]$D0_tpm,
#     WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64]$D4_tpm, method = "pearson")
# ggsave("/Users/KevinStein/Desktop/Deseq.gene.expression_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# ggplot(expression[sch9_D0_sum > 128 & sch9_D4_sum > 128], aes(sch9_D0_tpm, sch9_D4_tpm)) + geom_point() +
#   scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
# cor(expression[sch9_D0_sum > 128 & sch9_D4_sum > 128]$sch9_D0_tpm,
#     expression[sch9_D0_sum > 128 & sch9_D4_sum > 128]$sch9_D4_tpm, method = "pearson")
# 
# ggplot(expression[WT_D0_sum > 128 & sch9_D0_sum > 128], aes(WT_D0_tpm, sch9_D0_tpm)) + geom_point() +
#   scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
# cor(expression[WT_D0_sum > 128 & sch9_D0_sum > 128]$WT_D0_tpm,
#     expression[WT_D0_sum > 128 & sch9_D0_sum > 128]$sch9_D0_tpm, method = "pearson")
# 
# ggplot(expression[WT_D4_sum > 128 & sch9_D4_sum > 128], aes(WT_D4_tpm, sch9_D4_tpm)) + geom_point() +
#   scale_x_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   scale_y_log10(limits = c(1e-1,1e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
#   annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
# cor(expression[WT_D4_sum > 128 & sch9_D4_sum > 128]$WT_D4_tpm,
#     expression[WT_D4_sum > 128 & sch9_D4_sum > 128]$sch9_D4_tpm, method = "pearson")


# FigS2f: Gene ontology
GO <- read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/GO_translatome.csv", stringsAsFactors = T, header = T)
GO <- as.data.table(GO)
plot <- ggplot(GO[Organism == "yeast"], aes(x = reorder(Description, -foldchange), y = foldchange, fill = factor(Group))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("down","up"), values = c("#1F78B4", "#E31A1C"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_BP_translatomeYeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure 2 ###

# Fig 2a: Sequence motif
WT_stalling_peaks_ageDep <- readRDS("WT_stalling_peaks_ageDep.rds")
WT_peaks_all <- readRDS("WT_peaks_all.rds")
logomaker(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(WT_peaks_all_ageDep[WT_D4_pause > 6]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
#logomaker(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$motif3, type = "Logo", bg = bg, color_seed = 6)

# Fig 2b: 4 K/R
KR4of4_dt <- readRDS("doc/KR4of4_dt.rds")
plot <- ggplot(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale()) + coord_cartesian(ylim = c(0.5,2.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/WT_KR4of4.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Fig 2c: 6 K/R
KR6of6_dt <- readRDS("doc/KR6of6_dt.rds")
plot <- ggplot(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale()) + coord_cartesian(ylim = c(0.5,2.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/WT_KR6of6.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Fig 2d: Example
plot <- ggplot(data=sc_dtA[orf == "YOR272W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_rect(xmin = 242, xmax = 282, ymin = -Inf, ymax = Inf, fill = "#F2CF8C", alpha = 0.1) +
  geom_vline(xintercept = 262, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YOR272W_KR5of5_aa262.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YOR272W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T),
                                               S0 = movingAverage(((sch9_D0_1A_pause + sch9_D0_2A_pause) / 2), n=4, center=T),
                                               S4 = movingAverage(((sch9_D4_1A_pause + sch9_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 262, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  geom_line(aes(position, S0), color = "gray60", size = 1.25) +
  geom_line(aes(position, S4), color = "#00CED1", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(242,282))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YOR272W_YTM1_sch9.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



plot <- ggplot(data=sc_dtA[orf == "YOR272W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 262, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(242,282))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YOR272W_KR5of5_aa262_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# FigS3g: YEL056W example
plot <- ggplot(data=sc_dtA[orf == "YEL056W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 188, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expansion())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YEL056W_HAT2_Waa188.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(data=sc_dtA[orf == "YEL056W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T),
                                               S0 = movingAverage(((sch9_D0_1A_pause + sch9_D0_2A_pause) / 2), n=4, center=T),
                                               S4 = movingAverage(((sch9_D4_1A_pause + sch9_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 188, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  geom_line(aes(position, S0), color = "gray60", size = 1.25) +
  geom_line(aes(position, S4), color = "#00CED1", size = 1.25) +
  scale_x_continuous(expand = expansion(), limits = c(168,208))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YEL056W_HAT2_sch9.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Microscopy ###
microscopy <- as.data.table(read.csv("/Users/KevinStein/Desktop/Manuscript/Data/microscopy,gels/microscopy.csv", stringsAsFactors = T, header = T))
names(microscopy)
microscopy[, ID := paste0(Strain, "_", Reporter, "_", Day)]
microscopy[, strain_day := paste0(Strain, "_", Day)]
microscopy_GFPposRFPneg_fract <- as.data.table(summarySE(microscopy, measurevar="ratio", groupvars=c("Strain", "Day", "Reporter", "name")))
microscopy_GFPposRFPneg_fract[, strain_day := paste0(Strain, "_", Day)]

plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "R12" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "R12" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "ltn1", "ltn1", "rqc2", "rqc2", "hel2", "hel2")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_RQC_R12.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "K12" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "K12" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "ltn1", "ltn1", "rqc2", "rqc2", "hel2", "hel2")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_RQC_K12.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)


plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "GFP-RFP" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "GFP-RFP" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "ltn1", "ltn1", "rqc2", "rqc2", "hel2", "hel2")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_RQC_GFP-RFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)



plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "R12" & Strain %in% c("1WT", "sch9")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "R12" & Strain %in% c("1WT", "sch9")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "sch9", "sch9")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_sch9_R12.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "K12" & Strain %in% c("1WT", "sch9")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "K12" & Strain %in% c("1WT", "sch9")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "sch9", "sch9")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_sch9_K12.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "GFP-RFP" & Strain %in% c("1WT", "sch9")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "GFP-RFP" & Strain %in% c("1WT", "sch9")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "sch9", "sch9")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_sch9_GFP-RFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)


plot <- ggplot(microscopy_GFPposRFPneg_fract[Strain %in% c("hel2_GFP-RFP", "hel2_K12", "hel2_R12")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Strain %in% c("hel2_GFP-RFP", "hel2_K12", "hel2_R12")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  ylim(0,0.6)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.1,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_hel2_GFPRFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)


plot <- ggplot(microscopy_GFPposRFPneg_fract[Strain %in% c("ytm1_N", "ytm1_C")], 
               aes(x = name, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Strain %in% c("ytm1_N", "ytm1_C")], aes(name,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  ylim(0,0.6)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_ytm1_GFPRFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

plot <- ggplot(microscopy_GFPposRFPneg_fract[Strain %in% c("ytm1_N_ltn1", "ytm1_C_ltn1")], 
               aes(x = name, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Strain %in% c("ytm1_N_ltn1", "ytm1_C_ltn1")], aes(name,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.11_ytm1_ltn1_GFPRFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

plot <- ggplot(microscopy_GFPposRFPneg_fract[Strain %in% c("ytm1_N", "ytm1_C", "ytm1_N_ltn1", "ytm1_C_ltn1")], 
               aes(x = name, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Strain %in% c("ytm1_N", "ytm1_C", "ytm1_N_ltn1", "ytm1_C_ltn1")], aes(name,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(0.1, 1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.11_ytm1_ltn1ANDwt_GFPRFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)


plot <- ggplot(microscopy_GFPposRFPneg_fract[Strain %in% c("hat2_N")], 
               aes(x = name, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Strain %in% c("hat2_N")], aes(name,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  ylim(0,0.5)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/2021.06_hat2.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)




# plot <- ggplot(data=sc_dtA[orf == "YPL162C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=5, center=T), 
#                                                WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=5, center=T))]) + 
#   geom_rect(xmin = 53, xmax = 93, ymin = -Inf, ymax = Inf, fill = "#F2CF8C", alpha = 0.1) +
#   geom_vline(xintercept = 73, color = "gray75", linetype = "longdash", size = 1) +
#   geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
#   geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
#   scale_x_continuous(expand = expand_scale())
# plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
#   theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
# ggsave("/Users/KevinStein/Desktop/YPL162C_KR5of5_aa73.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
# 
# plot <- ggplot(data=sc_dtA[orf == "YPL162C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=5, center=T), 
#                                                WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=5, center=T))]) + 
#   geom_vline(xintercept = 73, color = "gray75", linetype = "longdash", size = 1) +
#   geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
#   geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
#   scale_x_continuous(expand = expand_scale(), limits = c(53,93))
# plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
#   theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
#         axis.text = element_text(size = 16, color = "black"))
# ggsave("/Users/KevinStein/Desktop/YPL162C_KR5of5_aa73_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure S3 ###

# FigS3a-c: tripeptide motifs
motif_dt <- readRDS("motif_dt.rds")
logomaker(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 6)
#logomaker(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 1.5]$motif, type = "Logo", bg = bg4_0, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio2_0 > 1 & D2_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 6)
#logomaker(motif_dt[count >= 100 & ratio2_0 > 1 & D2_pause > 1.5]$motif, type = "Logo", bg = bg4_0, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio4_2 > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_2, color_seed = 6)
#logomaker(motif_dt[count >= 100 & ratio4_2 > 1 & D4_pause > 1.5]$motif, type = "Logo", bg = bg4_2, color_seed = 6)

# FigS3d: Residue frequency
plot <- ggplot(aa_freq[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "A-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/A-siteFreq_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aa_freqP[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "P-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/P-siteFreq_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aa_freqE[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "E-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/E-siteFreq_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# FigS3e: Codon frequency
plot <- ggplot(codon_freq[codon != 'TAA' & codon != 'TAG' & codon != 'TGA'], aes(x = codon, y = log2(div), fill = residue)) + geom_col() +
  geom_hline(yintercept = 0, size = 0.2, color = 'black')
plot <- plot + theme_minimal(20, base_line_size = 0.5) + labs(y = "A-site frequency, log2", x = "Codon") +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 5, color = 'black'), axis.ticks.x = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 16, color = "black"), axis.ticks.y = element_line(color = "black"))
ggsave("/Users/KevinStein/Desktop/codonfreq_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

test <- matrix(c(aa_freq[residue == "C"]$stalls,
                 aa_freq[residue == "C"]$proteome,
                 (sum(aa_freq$stalls) - aa_freq[residue == "C"]$stalls),
                 (sum(aa_freq$proteome) - aa_freq[residue == "C"]$proteome)), nrow = 2)
fisher.test(test) # P = p-value < 2.2e-16; R: p-value = 0.006585; W: p-value = 0.02529; E: p-value = 1.549e-10
test <- matrix(c(aa_freqP[residue == "G"]$stalls,
                 aa_freqP[residue == "G"]$proteome,
                 (sum(aa_freqP$stalls) - aa_freqP[residue == "G"]$stalls),
                 (sum(aa_freqP$proteome) - aa_freqP[residue == "G"]$proteome)), nrow = 2)
fisher.test(test) # P = p-value = 0.01025; R: p-value = 3.356e-05; N: p-value = 1.559e-09; K: p-value = 1.033e-07; G: p-value = 3.447e-05
test <- matrix(c(aa_freqE[residue == "R"]$stalls,
                 aa_freqE[residue == "R"]$proteome,
                 (sum(aa_freqE$stalls) - aa_freqE[residue == "R"]$stalls),
                 (sum(aa_freqE$proteome) - aa_freqE[residue == "R"]$proteome)), nrow = 2)
fisher.test(test) # A = p-value = 0.02215; R: p-value = 0.05099; K: p-value = 0.00675

# FigS3f: Polybasic analysis

KR3of3_dt <- readRDS("doc/KR3of3_dt.rds")
KR5of5_dt <- readRDS("doc/KR5of5_dt.rds")

plot <- ggplot(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale()) + coord_cartesian(ylim = c(0.5,2.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/WT_KR3of3.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale()) + coord_cartesian(ylim = c(0.5,2.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/WT_KR5of5.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 3], 
               aes("3_0", ((WT_D0_1A_norm + WT_D0_2A_norm) / 2))) + 
  geom_boxplot(fill = "gray50", notch = T) +
  geom_boxplot(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 3], 
               aes("3_4", ((WT_D4_1A_norm + WT_D4_2A_norm) / 2)), fill = "#E7A427", notch = T) +
  geom_boxplot(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 4], 
               aes("4_0", ((WT_D0_1A_norm + WT_D0_2A_norm) / 2)), fill = "gray50", notch = T) +
  geom_boxplot(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 4], 
               aes("4_4", ((WT_D4_1A_norm + WT_D4_2A_norm) / 2)), fill = "#E7A427", notch = T) +
  geom_boxplot(data = KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 5], 
               aes("5_0", ((WT_D0_1A_norm + WT_D0_2A_norm) / 2)), fill = "gray50", notch = T) +
  geom_boxplot(data = KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 5], 
               aes("5_4", ((WT_D4_1A_norm + WT_D4_2A_norm) / 2)), fill = "#E7A427", notch = T) +
  geom_boxplot(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 6], 
               aes("6_0", ((WT_D0_1A_norm + WT_D0_2A_norm) / 2)), fill = "gray50", notch = T) +
  geom_boxplot(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 6], 
               aes("6_4", ((WT_D4_1A_norm + WT_D4_2A_norm) / 2)), fill = "#E7A427", notch = T) +
  coord_capped_cart(ylim = c(0,5), left = "both")
plot <- plot + theme_classic(18) + labs(y = "Norm. ribosome occupancy", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/polybasicStats_yeast.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)


KR3of3_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
KR3of3_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]
KR4of4_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
KR4of4_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]
KR5of5_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
KR5of5_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]
KR6of6_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
KR6of6_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]

wilcox.test(KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                        WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 3]$WT_D0_norm,
            KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                        WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 3]$WT_D4_norm, alternative = "t")
wilcox.test(KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                        WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 4]$WT_D0_norm,
            KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                        WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 4]$WT_D4_norm, alternative = "t")
wilcox.test(KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                        WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 5]$WT_D0_norm,
            KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                        WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 5]$WT_D4_norm, alternative = "t")
wilcox.test(KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                        WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 6]$WT_D0_norm,
            KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                        WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 25 & (polybasic - length) < -25 & adjusted > 0 & adjusted < 6]$WT_D4_norm, alternative = "t")

# plot <- ggplot(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
#                           WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
#                           polybasic > 25 & (polybasic - length) < -25,
#                         .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
#                           WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))]) + 
#   geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
#   stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
#                fun.args=list(conf.int=0.5), fill = 'gray30') +
#   stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
#                fun.args=list(conf.int=0.5), fill = '#FF7F00') +
#   stat_summary(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
#                                   WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
#                                   polybasic > 25 & (polybasic - length) < -25,
#                                 .(adjusted, WT0_6 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
#                                   WT4_6 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
#                aes(adjusted, WT0_6), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
#                fun.args=list(conf.int=0.5), fill = '#1F78B4') +
#   stat_summary(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
#                                   WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
#                                   polybasic > 25 & (polybasic - length) < -25,
#                                 .(adjusted, WT0_6 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
#                                   WT4_6 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
#                aes(adjusted, WT4_6), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
#                fun.args=list(conf.int=0.5), fill = '#E31A1C') +
#   stat_summary(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
#                                   WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
#                                   polybasic > 25 & (polybasic - length) < -25,
#                                 .(adjusted, WT0_6 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
#                                   WT4_6 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
#                aes(adjusted, WT4_6), fun.y = "mean", geom = "line", size=1.25, color = '#E31A1C') +
#   stat_summary(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
#                                   WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
#                                   polybasic > 25 & (polybasic - length) < -25,
#                                 .(adjusted, WT0_6 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
#                                   WT4_6 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
#                aes(adjusted, WT0_6), fun.y = "mean", geom = "line", size=1.25, color = '#1F78B4') +  
#   stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#FF7F00') +
#   stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray30') +
#   scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
# plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
#   theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
# ggsave("/Users/KevinStein/Desktop/WT_polybasic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Figure 3 ###

# Fig 3b: Lifespan and stalling
x <- mean(combined$foldchange)
y <- sd(combined$foldchange)

plot <- ggplot(combined[!is.na(Norm) & foldchange < x+y], aes("<1sd", log10(Norm_adj))) + 
  geom_hline(yintercept = 0.5456977, size = 1, color = "gray75", linetype = "dashed") +
  geom_violin(fill = "#E5F5E0") + geom_boxplot(fill = "#E5F5E0", width = 0.3) +
  geom_violin(data = combined[!is.na(Norm) & (foldchange > x+y) & (foldchange < x+3*y)], aes("1sd", log10(Norm_adj)), fill = "#A1D99B") +
  geom_boxplot(data = combined[!is.na(Norm) & (foldchange > x+y) & (foldchange < x+3*y)], aes("1sd", log10(Norm_adj)), fill = "#A1D99B", width = 0.3) +
  #geom_violin(data = combined[!is.na(Norm) & (foldchange > x+2*y) & (foldchange < x+3*y)], aes("2sd", log10(Norm_adj)), fill = "#41AB5D") +
  #geom_boxplot(data = combined[(foldchange > x+2*y) & (foldchange < x+3*y)], aes("2 SD", log10(Norm)), fill = "#41AB5D", width = 0.4) +
  geom_violin(data = combined[!is.na(Norm) & foldchange > x+3*y], aes("3sd", log10(Norm_adj)), fill = "#41AB5D") +
  geom_boxplot(data = combined[!is.na(Norm) & foldchange > x+3*y], aes("3sd", log10(Norm_adj)), fill = "#41AB5D", width = 0.3)
plot <- plot + theme_classic(20) + labs(y = "Norm. Lifespan", x = "GFP-R12 abundance") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Lifespan_Stalling_all.pdf", plot, width = 4, height = 6, dpi = 300, useDingbats = F)
wilcox.test(combined[!is.na(Norm) & foldchange < x+y]$Norm_adj, combined[!is.na(Norm) & (foldchange > x+y) & (foldchange < x+3*y)]$Norm_adj, alternative = 't')
wilcox.test(combined[!is.na(Norm) & foldchange < x+y]$Norm_adj, combined[!is.na(Norm) & (foldchange > x+3*y)]$Norm_adj, alternative = 't')
wilcox.test(combined[!is.na(Norm) & (foldchange > x+y) & (foldchange < x+3*y)]$Norm_adj, combined[!is.na(Norm) & (foldchange > x+3*y)]$Norm_adj, alternative = 't')

# plot <- ggplot(combined[!is.na(Norm) & foldchange > x+y], aes("2High", log10(Norm_adj))) + geom_violin(fill = "#238B45") + geom_boxplot(width = 0.4, fill = "#41AB5D") +
#   geom_violin(data = combined[!is.na(Norm) & foldchange < x+y], aes("1Low", log10(Norm_adj)), fill = "#A1D99B") +
#   geom_boxplot(data = combined[!is.na(Norm) & foldchange < x+y], aes("1Low", log10(Norm_adj)), width = 0.4, fill = "#C7E9C0")
# plot <- plot + theme_classic(20) + labs(y = "Norm. Lifespan", x = "GFP-R12 abundance") +
#   theme(axis.text = element_text(size = 20, color = "black"))
# ggsave("/Users/KevinStein/Desktop/Lifespan_Stalling.pdf", plot, width = 4, height = 6, dpi = 300, useDingbats = F)
# p <- wilcox.test(na.omit(combined[foldchange > x+y]$Norm_adj), na.omit(combined[foldchange < x+y]$Norm_adj), alternative = 't')
# p$p.value

# Fig 3c: WT vs sch9 age-dependent pause sites, metagene
WT_stalling_peaks_ageDep_dt <- readRDS("WT_stalling_peaks_ageDep_dt.rds")
plot <- ggplot(data = WT_stalling_peaks_ageDep_dt[WT_D0_1A_rpc_adjusted >= 1 & WT_D0_2A_rpc_adjusted >= 1 & 
                                                    WT_D4_1A_rpc_adjusted >= 1 & WT_D4_2A_rpc_adjusted >= 1 &
                                                    sch9_D0_1A_rpc_adjusted >= 1 & sch9_D0_2A_rpc_adjusted >= 1 & 
                                                    sch9_D4_1A_rpc_adjusted >= 1 & sch9_D4_2A_rpc_adjusted >= 1 &
                                                    peak > 25 & (peak - length) < -25,
                                                  .(adjusted, WT0 = movingAverage(WT_D0_norm, n=1, center=T),
                                                    WT4 = movingAverage(WT_D4_norm, n=1, center=T),
                                                    S0 = movingAverage(sch9_D0_norm, n=1, center=T),
                                                    S4 = movingAverage(sch9_D4_norm, n=1, center=T))]) +
  stat_summary(aes(adjusted, S0), fun.y = "mean", geom = "line", size=1.25, color = 'gray60') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') + 
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, S4), fun.y = "mean", geom = "line", size=1.25, color = '#00CED1') +
  scale_x_continuous(expand = expansion(), limits = c(-25,25))
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from pause site (codons)") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_sch9_ageDep.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Fig 3d: WT vs sch9 age-dependent pause sites, boxplot
WT_stalling_peaks <- readRDS("WT_stalling_peaks.rds")
plot <- ggplot(WT_stalling_peaks_ageDep[sch9_odds > 0], aes("1WT", WT_D4_pause)) + geom_boxplot(notch = T, fill = "#E7A427") +
  geom_boxplot(aes("2sch9", sch9_D4_pause), fill = "#00CED1", notch = T) + 
  coord_capped_cart(left = "both", ylim = c(0,5))
plot <- plot + theme_classic(20) + labs(y = "Pause score, Day 4", x = "") +
  theme(axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Day4_WTvsSch9dark.pdf", plot, width = 3, height = 6, dpi = 300, useDingbats = F)
test <- wilcox.test(WT_stalling_peaks_ageDep[sch9_D4_pause > 0]$WT_D4_pause, WT_stalling_peaks_ageDep[sch9_D4_pause > 0]$sch9_D4_pause, alternative = 't')
test$p.value

# Fig 3e: Logo at motifs with greater pausing in WT compared to sch9 (different motif_dt. See SequenceLogos.R)
logomaker(motif_dt[count >= 100 & ratio > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 6)

# Fig 3f: sch9 polybasic
plot <- ggplot(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  sch9_D0_1A_rpc_adjusted >= 0.1 & sch9_D0_2A_rpc_adjusted >= 0.1 & 
                                  sch9_D4_1A_rpc_adjusted >= 0.1 & sch9_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T),
                                  S0 = movingAverage(((sch9_D0_1A_norm + sch9_D0_2A_norm) / 2), n=3, center=T),
                                  S4 = movingAverage(((sch9_D4_1A_norm + sch9_D4_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, S0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray60') +
  stat_summary(aes(adjusted, S4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#00CED1') +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, S0), fun.y = "mean", geom = "line", size=1.25, color = 'gray60') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') + 
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, S4), fun.y = "mean", geom = "line", size=1.25, color = '#00CED1') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale()) + coord_cartesian(ylim = c(0.5,2.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/WT_sch9_KR4of4.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Fig 3g: sch9 polybasic
plot <- ggplot(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  sch9_D0_1A_rpc_adjusted >= 0.1 & sch9_D0_2A_rpc_adjusted >= 0.1 & 
                                  sch9_D4_1A_rpc_adjusted >= 0.1 & sch9_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T),
                                  S0 = movingAverage(((sch9_D0_1A_norm + sch9_D0_2A_norm) / 2), n=3, center=T),
                                  S4 = movingAverage(((sch9_D4_1A_norm + sch9_D4_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, S0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray60') +
  stat_summary(aes(adjusted, S4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#00CED1') +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, S0), fun.y = "mean", geom = "line", size=1.25, color = 'gray60') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') + 
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, S4), fun.y = "mean", geom = "line", size=1.25, color = '#00CED1') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale()) + coord_cartesian(ylim = c(0.5,2.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/WT_sch9_KR6of6.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# plot <- ggplot(WT_stalling_peaks[sch9_odds > 0], aes("1WT", WT_odds)) + geom_boxplot(notch = T, fill = "#E7A427") +
#   geom_boxplot(aes("2sch9", sch9_odds), fill = "#00CED1", notch = T) + 
#   coord_capped_cart(left = "both", ylim = c(0,5))
# plot <- plot + theme_classic(20) + labs(y = "odds ratio (Day 4 / Day 0)", x = "") +
#   theme(axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(),
#         axis.ticks.x = element_blank())
# ggsave("/Users/KevinStein/Desktop/Day4_WTvsSch9dark.pdf", plot, width = 3, height = 6, dpi = 300, useDingbats = F)
# test <- wilcox.test(WT_stalling_peaks[sch9_D4_pause > 0]$WT_odds, WT_stalling_peaks$sch9_odds, alternative = 't')
# test$p.value

### Metagene from start
# plot <- ggplot(data = sc_dtA[WT_D0_1A_rpc >= 1 & WT_D0_2A_rpc >= 1 & 
#                        WT_D2_1A_rpc >= 1 & WT_D2_2A_rpc >= 1 &
#                        WT_D4_1A_rpc >= 1 & WT_D4_2A_rpc >= 1 &
#                        WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
#                        WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
#                        WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & length >= 300]) + xlim(-7, 200) +
#   stat_summary(aes(position, WT_D0_pause, color = '#1F78B4'), fun.y = "mean", geom = "line", size= 1.25) +
#   stat_summary(aes(position, WT_D2_pause, color = '#999999'), fun.y = "mean", geom = "line", size= 1.25) +
#   stat_summary(aes(position, WT_D4_pause, color = '#E31A1C'), fun.y = "mean", geom = "line", size= 1.25) +
#   scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("#1F78B4", "#999999", "#E31A1C"), name = "")
# plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Codon position") +
#   theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
#         panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
# ggsave("/Users/KevinStein/Desktop/startMetagene_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

### Metagene from stop
# plot <- ggplot(data = sc_dtA[WT_D0_1A_rpc >= 1 & WT_D0_2A_rpc >= 1 & 
#                                WT_D2_1A_rpc >= 1 & WT_D2_2A_rpc >= 1 &
#                                WT_D4_1A_rpc >= 1 & WT_D4_2A_rpc >= 1 &
#                                WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
#                                WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
#                                WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & length >= 300,
#                              .(stopdist, WT0 = WT_D0_pause, WT2 = WT_D2_pause,
#                                WT4 = WT_D4_pause)]) + xlim(-50, 7) +
#   stat_summary(aes(stopdist, WT4, color = '#E31A1C'), fun.y = "mean", geom = "line", size= 1.25) +
#   stat_summary(aes(stopdist, WT2, color = '#999999'), fun.y = "mean", geom = "line", size= 1.25) +
#   stat_summary(aes(stopdist, WT0, color = '#1F78B4'), fun.y = "mean", geom = "line", size= 1.25) +
#   scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("#1F78B4", "#999999", "#E31A1C"), name = "")
# plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Codon position") +
#   theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
#         panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
# ggsave("/Users/KevinStein/Desktop/stopMetagene_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


temp1 <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Yeast/YeastProteomeCategorization/Localization_Uniprot.csv", header = T)
temp1 <- as.data.table(temp1)
i <- cbind(match(WT_stalling_peaks_ageDep$orf, temp1$orf))
WT_stalling_peaks_ageDep <- cbind(WT_stalling_peaks_ageDep, name = temp1[i]$Names)
write.csv(WT_stalling_peaks_ageDep, "/Users/KevinStein/Desktop/yeast.csv")
