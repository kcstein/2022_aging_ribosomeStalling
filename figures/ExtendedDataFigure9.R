### Extended Data Figure 9 ###


# Fig S9a: Residue frequency
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


# Fig S9b: Codon frequency
plot <- ggplot(codon_freq[codon != 'TAA' & codon != 'TAG' & codon != 'TGA'], aes(x = codon, y = log2(div), fill = residue)) + geom_col() +
  geom_hline(yintercept = 0, size = 0.2, color = 'black')
plot <- plot + theme_minimal(20, base_line_size = 0.5) + labs(y = "A-site frequency, log2", x = "Codon") +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 5, color = 'black'), axis.ticks.x = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 16, color = "black"), axis.ticks.y = element_line(color = "black"))
ggsave("/Users/KevinStein/Desktop/codonfreq_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S9c-e: Sequence logos
motif_dt <- readRDS("motif_dt.rds")
logomaker(motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_6_1 > 1 & D6_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_12_6 > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg_noW, color_seed = 6)


# Fig S9f: K12C11.6 example, W in P-site
plot <- ggplot(data=N2_dtA[orf == "K12C11.6", .(position, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=3, center=T), 
                                               D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=3, center=T))]) + 
  #geom_vline(xintercept = 62, color = "gray50", linetype = "longdash", size = 1) +
  geom_line(aes(position, D12), color = "#E7298A", size = 1.25) +
  geom_line(aes(position, D1), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expansion())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/K12C11.6_Waa62.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S9g: Polybasic analysis
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

