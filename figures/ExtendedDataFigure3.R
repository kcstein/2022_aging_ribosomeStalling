### Extended Data Figure 3 ###


# Fig S3a: Residue pausing in yeast treated with 3AT
plot <- ggplot(pauseMean, aes(WT, AT, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5)
plot <- plot + theme_classic(20) + labs(y = "AT pause score", x = "WT pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Green_ATvsWT.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S3b: Residue frequency in A-site at pause sites in yeast treated with 3AT
plot <- ggplot(aa_freq[residue != 'X'], aes(residue, log2(div), fill = residue)) + geom_col() + 
  geom_hline(yintercept = 0, color = 'black', size = 0.2) + coord_cartesian(ylim = c(-4,2.5))
plot <- plot + theme_classic(20) + labs(y = "A-site frequency, log2", x = "Residue") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/A-siteFreq_Green.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S3c: Yeast GO of genes with pausing
GO <- read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/GO_stalling.csv", stringsAsFactors = T, header = T)
GO <- as.data.table(GO)
plot <- ggplot(GO[Organism == "yeast"], aes(x = reorder(Description, Fold.Enrichment), y = Fold.Enrichment, fill = factor(Category))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("BP","CC"), values = c("#A6CEE3", "#B2DF8A"), name = "")
plot <- plot + theme_classic(12) + labs(y = "Fold enrichment", x = "") +
  theme(legend.position = c(0.8,0.2), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_stallingYeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S3d: Worm GO of genes with pausing
plot <- ggplot(GO[Organism == "worms"], aes(x = reorder(Description, Fold.Enrichment), y = Fold.Enrichment, fill = factor(Category))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("BP","CC"), values = c("#A6CEE3", "#B2DF8A"), name = "")
plot <- plot + theme_classic(12) + labs(y = "Fold enrichment", x = "") +
  theme(legend.position = c(0.8,0.2), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_stallingWorm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S3e: Association of stalling and co-translational ubiquitination
temp <- data.table(Category = c("No", "Stall"), 
                   FractionUbiq = c(length(ubiq_subset[!ubiq_subset$YORF %in% stall$orf]$YORF) /
                                      length(proteome_orfs[!proteome_orfs$orf %in% stall$orf]$orf),
                                    length(ubiq_subset[ubiq_subset$YORF %in% stall$orf]$YORF) /
                                      length(stall$orf)))
plot <- ggplot(temp, aes(x = Category, y=FractionUbiq, fill = Category)) + geom_col(color = "black", size = 0.5) + scale_y_continuous(limits = c(0,0.15), expand = c(0.0008,0.0008)) +
  scale_fill_manual(limits = c("No","Stall"), values = c("#80B1D3", "#FDB462"), name = "")
plot <- plot + theme_classic(18) + labs(y = "Fraction of ubiquitinated\nproteins in dataset", x = "") +
  theme(legend.position = "none", legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Stalling_DuttlerUbiq_0.002388.pdf", plot, width = 3.2, height = 4, dpi = 300, useDingbats = F)


# Fig S3f: Codon frequency
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


# Fig S3g-i: tripeptide motifs
motif_dt <- readRDS("motif_dt.rds")
logomaker(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio2_0 > 1 & D2_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio4_2 > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_2, color_seed = 6)


