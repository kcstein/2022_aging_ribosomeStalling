### Figure 2 ###

# Fig 2a: Sequence motif
WT_stalling_peaks_ageDep <- readRDS("WT_stalling_peaks_ageDep.rds")
logomaker(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$motif3, type = "EDLogo", bg = bg, color_seed = 6)


# Fig 2b: Residue frequency
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


# Fig 2c: 4 K/R
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


# Fig 2d: 6 K/R
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


# Fig 2e: YTM1
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
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 262, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(242,282))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YOR272W_KR5of5_aa262_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig 2f: YTM1 disome
plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YOR272W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=4, center=T),
                                                                            occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=4, center=T))]) + 
  geom_rect(xmin = 242, xmax = 282, ymin = -Inf, ymax = Inf, fill = "#AFEEEE", alpha = 0.01) +
  geom_vline(xintercept = 262, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  geom_line(color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/Ytm1_disomeMonosome.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YOR272W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=3, center=T),
                                                                            occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=3, center=T))]) + 
  geom_vline(xintercept = 262, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  geom_line(color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(242,282)) + coord_cartesian(ylim = c(0,1.5))
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/Ytm1_disomeMonosome_inset.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

