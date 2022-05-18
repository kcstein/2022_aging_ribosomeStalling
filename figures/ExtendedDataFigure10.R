### Extended Data Figure 10 ###

# Fig S10a: Association of pausing and aggregation
temp <- data.table(Category = c("No", "Stall"), 
                   FractionAggregation = c(((length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                                               length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)) /
                                              (length(hartl_proteome_orfs$WBgene) - 4)),
                                           (length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)) /
                                             length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene)))
plot <- ggplot(temp, aes(x = Category, y=FractionAggregation, fill = Category)) + geom_col(color = "black", size = 0.5) + scale_y_continuous(limits = c(0,0.15), expand = c(0.0005,0.0005)) +
  scale_fill_manual(limits = c("No","Stall"), values = c("#80B1D3", "#FDB462"), name = "")
plot <- plot + theme_classic(18) + labs(y = "Fraction of aggregated\nproteins in dataset", x = "") +
  theme(legend.position = "none", legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Stalling_kenyonAggregation_p5.6e-22.pdf", plot, width = 3.2, height = 4, dpi = 300, useDingbats = F)


# Fig S10b: Enrichment of basic residues
plot <- ggplot(agg_polybasic_dt[datasource == "Hartl"], aes(x = residues, y = odds, fill = stall)) + 
  geom_col(color = "black", size = 0.5, position = "dodge") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(limits = c("a","y"), labels = c("all", "yes"), values = c("gray50", "#00CED1"), name = "Pause site") +
  scale_y_continuous(expand = c(0.0008,0.0008)) +
  scale_x_discrete(labels = c("4 Lys/Arg", "2 Arg"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio (relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.1,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_HartlAggregation.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)

plot <- ggplot(agg_polybasic_dt[datasource == "Kenyon"], aes(x = residues, y = odds, fill = stall)) + 
  geom_col(color = "black", size = 0.5, position = "dodge") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(limits = c("a","y"), labels = c("all", "yes"), values = c("gray50", "#00CED1"), name = "Pause site") +
  scale_y_continuous(expand = c(0.0008,0.0008)) +
  scale_x_discrete(labels = c("4 Lys/Arg", "2 Arg"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio (relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.1,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_KenyonAggregation.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)


# Fig S10c: GO of genes with pausing and aggregation
plot <- ggplot(GO_Kenyon_subset, aes(x = reorder(Description, div), y = Fold.Enrichment, fill = factor(Group))) + 
  geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("all","stall"), values = c("#FB9A99", "#E31A1C"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = c(0.85,0.15), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_Kenyon_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S10d: Aggregation of tRNA synthetases
tRNA_aggregation <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/tRNA_stall.csv", header = TRUE, stringsAsFactors = TRUE))
plot <- ggplot(tRNA_aggregation1, aes(age,log2.ratio.lh, color = gene, group = gene)) + geom_line(size = 1) + geom_point(size = 2)
plot <- plot + theme_classic(18) + labs(y = "Aggregate abundance, log2", x = "Age (days)") +
  theme(legend.position = c(0.75,0.4), legend.text = element_text(size = 12), legend.title = element_blank(), legend.background = element_blank(),
        axis.text = element_text(color = "black", size = 14),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/tRNA_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S10e: RQC ortholog expression
wormRQCexpression <- as.data.table(read.csv("/Users/KevinStein/Desktop/N2_RQCexpression.csv", header = TRUE))
plot <- ggplot(wormRQCexpression, aes(orf, log2FoldChange)) + geom_col() +
  labs(y = "Fold change (Day12 / Day1), log2", x = "")
plot <- plot + theme_classic(15) +
  theme(axis.line.x = element_blank(), axis.text.y = element_text(size = 12, color = "black"), axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 12, color = 'black'))
ggsave("/Users/KevinStein/Desktop/RQCexpression_worm.pdf", plot, width = 4, height = 6, dpi = 300, useDingbats = F)


# Fig S10f: Aggregation of RQC orthologs
rqc_aggregation <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/rqc_aggregation.csv", header = TRUE, stringsAsFactors = TRUE))
plot <- ggplot(rqc_aggregation1[gene != "C52E12.1"], aes(age,log2.ratio.lh, color = gene, group = gene)) + geom_line(size = 1) + geom_point(size = 2)
plot <- plot + theme_classic(18) + labs(y = "Aggregate abundance, log2", x = "Age (days)", color = "") +
  theme(legend.position = c(0.25,0.85), legend.text = element_text(size = 12), legend.title = element_blank(), legend.background = element_blank(),
        axis.text = element_text(color = "black", size = 14),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/rqc_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

