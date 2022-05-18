library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')

# yeast data
sc_dtA <- readRDS("doc/sc_dtA.rds") 
sc_fishers <- readRDS("doc/sc_fishers.rds")
WT_peaks_all <- readRDS("doc/WT_peaks_all.rds")
sch9_peaks_all <- readRDS("sch9_peaks_all.rds")
WT_stalling_peaks_ageDep_dt <- readRDS("WT_stalling_peaks_ageDep_dt.rds")

# worm data
N2_dtA <- readRDS("doc/N2_dtA.rds") 
N2_fishers <- readRDS("doc/N2_fishers.rds")
stalling_peaks_all <- readRDS("stalling_peaks_all.rds")
stalling_peaks_ageDep_dt <- readRDS("stalling_peaks_ageDep_dt.rds")


### Figure 1

# Fig 1c: Global stalling
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



# Fig 1d: Volcano of age-dependent stalls
plot <- ggplot(N2_fishers[odds < Inf & odds > 0 & padj != 0]) + 
  geom_point(aes(log2(odds), -log10(padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = stalling_peaks_all[D12_pause > N2_D6_2A_pause & padj != 0], aes(log2(odds), -log10(padj), color = "#E7298A"), alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous() + 
  scale_color_manual(labels = c(""), values = c("#E7298A"), name = "")
  #scale_color_manual(labels = c("Pausing","Age-dep. pausing"), values = c("#9E9AC8", "#FF7F00"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D12/D1)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_worms.tiff", plot, width = 6, height = 4, dpi = 300)

plot <- ggplot(sc_fishers[WT_odds < Inf & WT_odds > 0 & WT_padj != 0]) + 
  geom_point(aes(log2(WT_odds), -log10(WT_padj)), color = "gray80", alpha = 0.5, size = 2) +
  geom_point(data = WT_peaks_all[WT_D4_pause > WT_D2_pause & WT_padj != 0], aes(log2(WT_odds), -log10(WT_padj), color = "#E7A427"), alpha = 0.6, size = 2) +
  scale_x_continuous() + scale_y_continuous(limits = c(0,200)) + 
  scale_color_manual(labels = c(""), values = c("#E7A427"), name = "")
  #scale_color_manual(labels = c("Pausing","Age-dep. pausing"), values = c("#9E9AC8", "#FF7F00"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Relative pausing (log2 odds ratio D4/D0)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Volcano_stalling_yeast.tiff", plot, width = 6, height = 4, dpi = 300)


# Fig 1e: Metagene of age-dependent pausing sites
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

