### Extended Data Figure 8 ###


# Fig S8a: RQC flux in various categories of yeast strains
plot <- ggplot(data = combined[Norm_adj > x1+y1], aes("2CLS", log2(foldchange))) + geom_boxplot(fill = "gray60", notch = T) +
  geom_boxplot(data = combined[foldchange > x+3*y], aes("1RQC", log2(foldchange)), fill = "#B2DF8A", notch = T) +
  geom_boxplot(data = stall[stall$ORF %in% tor$gene], aes("3TOR", log2(foldchange)), fill = "gray60", notch = T) +
  geom_boxplot(data = stall[stall$ORF %in% RLS$gene], aes("5RLS", log2(foldchange)), fill = "gray60", notch = T) +
  geom_boxplot(data = combined[(!combined$ORF %in% tor$gene) & Norm_adj > x1+y1], aes("4CLS_noTOR", log2(foldchange)), fill = "gray60", notch = T) +
  geom_point(data = stall[ORF == "YHR205W"], aes("2CLS", log2(foldchange)), color = '#E7298A', size = 2) +
  labs(y = "RQC flux (GFP-R12 fold change, log2)", x = "")
plot <- plot + theme_classic(15) + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/RQCflux.pdf", plot, width = 5, height = 5, dpi = 300, useDingbats = F)


# Fig S8b: Correlation of datasets
library(corrplot)
library(tidyr)
names(expression)
M <- cor(expression[, c(35,37,39,41,43,45,23,25,27,29,31,33)])
head(round(M, 2))
corrplot(M, method="color",  
         type="upper", 
         tl.col="black", tl.srt=45)
pdf(file = "/Users/KevinStein/Desktop/WT_sch9_correlation.pdf", width = 6.5, height = 5.5, useDingbats = F)
corrplot(M, method="color",  
         type="upper", 
         tl.col="black", tl.srt=45, cl.align.text = 'l', addCoef.col = "white", number.cex = 0.7
         )
dev.off()


# Fig S8c: WT vs sch9 age-dependent pause sites, metagene
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

plot <- ggplot(WT_stalling_peaks_ageDep[sch9_odds > 0], aes("1WT", WT_D4_pause)) + geom_boxplot(notch = T, fill = "#E7A427") +
  geom_boxplot(aes("2sch9", sch9_D4_pause), fill = "#00CED1", notch = T) + 
  coord_capped_cart(left = "both", ylim = c(0,5))
plot <- plot + theme_classic(20) + labs(y = "Pause score, Day 4", x = "") +
  theme(axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Day4_WTvsSch9dark.pdf", plot, width = 3, height = 6, dpi = 300, useDingbats = F)
test <- wilcox.test(WT_stalling_peaks_ageDep[sch9_D4_pause > 0]$WT_D4_pause, WT_stalling_peaks_ageDep[sch9_D4_pause > 0]$sch9_D4_pause, alternative = 't')
test$p.value


# Fig S8d: Logo at motifs with greater pausing in WT compared to sch9 (different motif_dt. See SequenceLogos.R)
logomaker(motif_dt[count >= 100 & ratio > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 6)


# Fig S8e: sch9 polybasic KR6
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


# Fig S8f: YTM1 and HAT2 inset with sch9 KO data
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


# Fig S8i: RQC expression
RQCexpression <- as.data.table(read.csv("/Users/KevinStein/Desktop/yeast_RQCexpression.csv", header = TRUE))

plot <- ggplot(RQCexpression, aes(gene, log2FoldChange, fill = ID)) + geom_col(position = "dodge") +
  labs(y = "Fold change (Day4 / Day0), log2", x = "", fill = "") + 
  scale_fill_discrete(limits = c("a", "b"), labels = c("WT", "sch9"))
plot <- plot + theme_classic(15) +
  theme(legend.position = c(.8,.3), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=12),
        axis.line.x = element_blank(), axis.text = element_text(size = 15, color = "black"), axis.ticks.x = element_blank())
ggsave("/Users/KevinStein/Desktop/RQCexpression_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

