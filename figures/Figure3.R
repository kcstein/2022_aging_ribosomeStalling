### Figure 3 ###


# Fig 3b: R12 microscopy
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
ggsave("/Users/KevinStein/Desktop/RQC_R12.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)


# Fig 3c: Lifespan and stalling
x <- mean(combined$foldchange)
y <- sd(combined$foldchange)

plot <- ggplot(combined[!is.na(Norm) & foldchange < x+y], aes("<1sd", log10(Norm_adj))) + 
  geom_hline(yintercept = 0.5456977, size = 1, color = "gray75", linetype = "dashed") +
  geom_violin(fill = "#E5F5E0") + geom_boxplot(fill = "#E5F5E0", width = 0.3) +
  geom_violin(data = combined[!is.na(Norm) & (foldchange > x+y) & (foldchange < x+3*y)], aes("1sd", log10(Norm_adj)), fill = "#A1D99B") +
  geom_boxplot(data = combined[!is.na(Norm) & (foldchange > x+y) & (foldchange < x+3*y)], aes("1sd", log10(Norm_adj)), fill = "#A1D99B", width = 0.3) +
  geom_violin(data = combined[!is.na(Norm) & foldchange > x+3*y], aes("3sd", log10(Norm_adj)), fill = "#41AB5D") +
  geom_boxplot(data = combined[!is.na(Norm) & foldchange > x+3*y], aes("3sd", log10(Norm_adj)), fill = "#41AB5D", width = 0.3)
plot <- plot + theme_classic(20) + labs(y = "Norm. Lifespan", x = "GFP-R12 abundance") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Lifespan_Stalling_all.pdf", plot, width = 4, height = 6, dpi = 300, useDingbats = F)
wilcox.test(combined[!is.na(Norm) & foldchange < x+y]$Norm_adj, combined[!is.na(Norm) & (foldchange > x+y) & (foldchange < x+3*y)]$Norm_adj, alternative = 't')
wilcox.test(combined[!is.na(Norm) & foldchange < x+y]$Norm_adj, combined[!is.na(Norm) & (foldchange > x+3*y)]$Norm_adj, alternative = 't')
wilcox.test(combined[!is.na(Norm) & (foldchange > x+y) & (foldchange < x+3*y)]$Norm_adj, combined[!is.na(Norm) & (foldchange > x+3*y)]$Norm_adj, alternative = 't')


# Fig 3d: sch9 polybasic KR4
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


# Fig 3e: sch9 microscopy
plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "R12" & Strain %in% c("1WT", "sch9")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "R12" & Strain %in% c("1WT", "sch9")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "sch9", "sch9")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/sch9_R12.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "K12" & Strain %in% c("1WT", "sch9")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "K12" & Strain %in% c("1WT", "sch9")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "sch9", "sch9")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/sch9_K12.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "GFP-RFP" & Strain %in% c("1WT", "sch9")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "GFP-RFP" & Strain %in% c("1WT", "sch9")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "sch9", "sch9")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/sch9_GFP-RFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

