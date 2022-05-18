### Extended Data Figure 4 ###


# Fig S4a: Residue enrichment around disome positions
logomaker(motif_dt[count >= 5 & ratio > 1]$motif, type = "EDLogo", bg = bg, color_seed = 6)


# Fig S4b: HAT2
plot <- ggplot(data=sc_dtA[orf == "YEL056W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 188, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expansion())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YEL056W_HAT2_Waa188.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S4c: HAT2 disome
plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YEL056W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=4, center=T),
                                                                            occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=4, center=T))]) + 
  geom_rect(xmin = 168, xmax = 208, ymin = -Inf, ymax = Inf, fill = "#AFEEEE", alpha = 0.01) +
  geom_vline(xintercept = 188, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  geom_line(color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/Hat2_disomeMonosome.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S4d: Maximal disome positions
plot <- ggplot(data = disome_peaks_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                                WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 &
                                disome > 25 & (disome - length) < -25,
                              .(adjusted, WT0 = movingAverage(WT_D0_norm, n=1, center=T),
                                WT2 = movingAverage(WT_D2_norm, n=1, center=T),
                                WT4 = movingAverage(WT_D4_norm, n=1, center=T))]) +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray40", "#F2CF8C", "#E7A427"), name = "") +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from disome site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_disomepositions_Age.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = disome_peaks_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                                        WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 &
                                        disome > 25 & (disome - length) < -25,
                                      .(adjusted, WT0 = movingAverage(WT_D0_norm, n=1, center=T),
                                        WT2 = movingAverage(WT_D2_norm, n=1, center=T),
                                        WT4 = movingAverage(WT_D4_norm, n=1, center=T))]) +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  #stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = '#F2CF8C') +  
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray40", "#F2CF8C", "#E7A427"), name = "") +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.3))
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from disome site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_disomepositions_Age_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S4e: Inhibitory codons
plot <- ggplot(data = dicodon_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                                  WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                                  codon_position > 25 & (codon_position - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/dicodon.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = dicodon_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                           WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                           codon_position > 25 & (codon_position - length) < -25,
                         .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                           WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/dicodon_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S4f: RKK motif
plot <- ggplot(data = RKK_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/RKK.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = RKK_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.2))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/RKK_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S4g: GGG motif
plot <- ggplot(data = GGG_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/GGG.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = GGG_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.2))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/GGG_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S4h: PPG motif
plot <- ggplot(data = PPG_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/PPG.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = PPG_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.2))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/PPG_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

