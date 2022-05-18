### Aggregation of tRNA synthetases
tRNA <- c("WBGene00001094", "WBGene00006617", "WBGene00006945", 
          "WBGene00006936", "WBGene00004679", "WBGene00003815", 
          "WBGene00003415", "WBGene00003073", "WBGene00001336")
tRNA_aggregation <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/tRNA_stall.csv", header = TRUE, stringsAsFactors = TRUE))

names(tRNA_aggregation)
tRNA_aggregation1 <- gather(tRNA_aggregation, age, log2.ratio.lh, "Day1":"Day17", factor_key=TRUE)
tRNA_aggregation1 <- as.data.table(tRNA_aggregation1)
setkeyv(tRNA_aggregation1, c("uniprot"))

plot <- ggplot(tRNA_aggregation1, aes(age,log2.ratio.lh, color = gene, group = gene)) + geom_line(size = 1) + geom_point(size = 2)
plot <- plot + theme_classic(18) + labs(y = "Aggregate abundance, log2", x = "Age (days)") +
  theme(legend.position = c(0.75,0.4), legend.text = element_text(size = 12), legend.title = element_blank(), legend.background = element_blank(),
        axis.text = element_text(color = "black", size = 14),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/tRNA_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Aggregation of RQC factors
worm_rqc <- c("WBGene00010556", "WBGene00010620", "WBGene00022350", 
          "WBGene00016888", "WBGene00021831")
rqc_aggregation <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/rqc_aggregation.csv", header = TRUE, stringsAsFactors = TRUE))

names(rqc_aggregation)
rqc_aggregation1 <- gather(rqc_aggregation, age, log2.ratio.lh, "Day1":"Day17", factor_key=TRUE)
rqc_aggregation1 <- as.data.table(rqc_aggregation1)
setkeyv(rqc_aggregation1, c("uniprot"))

plot <- ggplot(rqc_aggregation1[gene != "C52E12.1"], aes(age,log2.ratio.lh, color = gene, group = gene)) + geom_line(size = 1) + geom_point(size = 2)
plot <- plot + theme_classic(18) + labs(y = "Aggregate abundance, log2", x = "Age (days)", color = "") +
  theme(legend.position = c(0.25,0.85), legend.text = element_text(size = 12), legend.title = element_blank(), legend.background = element_blank(),
        axis.text = element_text(color = "black", size = 14),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/rqc_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)