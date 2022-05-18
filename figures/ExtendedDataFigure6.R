### Extended Data Figure 6 ###


# Fig S6f: K12 and GFP-RFP microscopy
microscopy <- as.data.table(read.csv("/Users/KevinStein/Desktop/Manuscript/Data/microscopy,gels/microscopy.csv", stringsAsFactors = T, header = T))
names(microscopy)
microscopy[, ID := paste0(Strain, "_", Reporter, "_", Day)]
microscopy[, strain_day := paste0(Strain, "_", Day)]
microscopy_GFPposRFPneg_fract <- as.data.table(summarySE(microscopy, measurevar="ratio", groupvars=c("Strain", "Day", "Reporter", "name")))
microscopy_GFPposRFPneg_fract[, strain_day := paste0(Strain, "_", Day)]


plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "K12" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "K12" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "ltn1", "ltn1", "rqc2", "rqc2", "hel2", "hel2")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/RQC_K12.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

plot <- ggplot(microscopy_GFPposRFPneg_fract[Reporter == "GFP-RFP" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], 
               aes(x = strain_day, y = ratio, fill = as.factor(Day))) + geom_bar(stat="identity", size = 0.5, color = "black", position = "dodge") +
  geom_point(data = microscopy[Reporter == "GFP-RFP" & Strain %in% c("1WT", "2ltn1", "3rqc2", "4hel2")], aes(strain_day,ratio), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2, position = "dodge") + 
  scale_fill_manual(limits = c("0", "4"), labels = c("Day 0", "Day 4"), values = c("gray40", "#E7A427"), name = "") + 
  scale_x_discrete(labels = c("WT", "WT", "ltn1", "ltn1", "rqc2", "rqc2", "hel2", "hel2")) + ylim(0,0.8)
plot <- plot + theme_classic(20) + labs(y = "Fraction of cells", x = "") + 
  theme(legend.position = c(.7,1.1), legend.justification = c("left", "top"), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 16))
ggsave("/Users/KevinStein/Desktop/RQC_GFP-RFP.pdf", plot, width = 4.5, height = 6, dpi = 300, useDingbats = F)

