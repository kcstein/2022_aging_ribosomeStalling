### Extended Data Figure 7 ###


# Fig S7b: YTM1 microscopy
microscopy <- as.data.table(read.csv("/Users/KevinStein/Desktop/Manuscript/Data/microscopy,gels/microscopy.csv", stringsAsFactors = T, header = T))
names(microscopy)
microscopy[, ID := paste0(Strain, "_", Reporter, "_", Day)]
microscopy[, strain_day := paste0(Strain, "_", Day)]
microscopy_GFPposRFPneg_fract <- as.data.table(summarySE(microscopy, measurevar="ratio", groupvars=c("Strain", "Day", "Reporter", "name")))
microscopy_GFPposRFPneg_fract[, strain_day := paste0(Strain, "_", Day)]

plot <- ggplot(microscopy_GFPposRFPneg_fract[Strain %in% c("ytm1_N", "ytm1_C", "ytm1_N_ltn1", "ytm1_C_ltn1")], 
               aes(x = name, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Strain %in% c("ytm1_N", "ytm1_C", "ytm1_N_ltn1", "ytm1_C_ltn1")], aes(name,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(0.1, 1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/ytm1_ltn1ANDwt_GFPRFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)


# Fig S7d: Hel2 GFP+/RFP+ microscopy
plot <- ggplot(microscopy_GFPposRFPneg_fract[Strain %in% c("hel2_GFP-RFP", "hel2_K12", "hel2_R12")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Strain %in% c("hel2_GFP-RFP", "hel2_K12", "hel2_R12")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  ylim(0,0.6)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.1,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/hel2_GFPRFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)


# Fig S7e: Stalling reporter
plot <- ggplot(aes(position, occupancy), data=GHA_dt[orf == "GFP_AAA12_HIS3", .(position, occupancy = movingAverage(GHA0_pause, n=5, center=T),
                                                                       occupancy2 = movingAverage(GHA4_pause, n=5, center=T))]) + 
  geom_line(color = "gray40", size = 1.25) +
  geom_line(aes(position, occupancy2), color = '#E7A427', size = 1.25)
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (pause score)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/K12reporter_pause.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=GH_dt[orf == "YJR009C", .(position, occupancy = movingAverage(GH0_pause, n=1, center=T),
                                                                    occupancy2 = movingAverage(GH4_pause, n=1, center=T))]) + 
  geom_vline(xintercept = 310, color = "gray50", linetype = "dashed") +
  geom_line(aes(position, occupancy2), color = 'red', size = 1.25) +
  geom_line(color = "gray30", size = 1.25)  
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (pause score)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/GHreporter_pause.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

