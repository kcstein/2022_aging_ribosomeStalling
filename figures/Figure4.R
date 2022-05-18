### Figure 4 ###

# Fig 4a: Sequence logo
stalling_peaks_ageDep <- readRDS("stalling_peaks_ageDep.rds")
logomaker(stalling_peaks_ageDep[D12_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
#logomaker(stalling_peaks_all_ageDep[D12_pause > 10]$motif3, type = "EDLogo", bg = bg, color_seed = 6)


# Fig 4b: A-site residue frequency
plot <- ggplot(aa_freq[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-0.9,0.6))
plot <- plot + theme_classic(20) + labs(y = "A-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/A-siteFreq_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

test <- matrix(c(aa_freq[residue == "R"]$stalls,
                 aa_freq[residue == "R"]$proteome,
                 (sum(aa_freq$stalls) - aa_freq[residue == "R"]$stalls),
                 (sum(aa_freq$proteome) - aa_freq[residue == "R"]$proteome)), nrow = 2)
fisher.test(test) # R: p-value = 0.004836


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

