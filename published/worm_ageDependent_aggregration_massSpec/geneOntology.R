### Gene ontology
GO_agg <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/GO/GO_summary_worm_aggregation.csv", header = TRUE, stringsAsFactors = TRUE))
names(GO_agg)


GO_Hartl <- GO_agg[Dataset == "Hartl" & Group == "all"]
GO_Hartl_stall <- GO_agg[Dataset == "Hartl" & Group == "stall"]
GO_Kenyon <- GO_agg[Dataset == "Kenyon" & Group == "all"]
GO_Kenyon_stall <- GO_agg[Dataset == "Kenyon" & Group == "stall"]

GO_Hartl_common <- GO_Hartl[GO_Hartl$Term %in% GO_Hartl_stall$Term]
GO_Kenyon_common <- GO_Kenyon[GO_Kenyon$Term %in% GO_Kenyon_stall$Term]

i <- cbind(match(GO_Hartl_common$Term, GO_Hartl_stall$Term))
GO_Hartl_common <- cbind(GO_Hartl_common, foldEnrich_stall = GO_Hartl_stall[i]$Fold.Enrichment)
GO_Hartl_common <- cbind(GO_Hartl_common, count_stall = GO_Hartl_stall[i]$Count)
GO_Hartl_common <- cbind(GO_Hartl_common, genes_stall = GO_Hartl_stall[i]$Genes)
GO_Hartl_common <- cbind(GO_Hartl_common, benjamini_stall = GO_Hartl_stall[i]$Benjamini)
GO_Hartl_common <- cbind(GO_Hartl_common, pvalue_stall = GO_Hartl_stall[i]$PValue)

i <- cbind(match(GO_Kenyon_common$Term, GO_Kenyon_stall$Term))
GO_Kenyon_common <- cbind(GO_Kenyon_common, foldEnrich_stall = GO_Kenyon_stall[i]$Fold.Enrichment)
GO_Kenyon_common <- cbind(GO_Kenyon_common, count_stall = GO_Kenyon_stall[i]$Count)
GO_Kenyon_common <- cbind(GO_Kenyon_common, genes_stall = GO_Kenyon_stall[i]$Genes)
GO_Kenyon_common <- cbind(GO_Kenyon_common, benjamini_stall = GO_Kenyon_stall[i]$Benjamini)
GO_Kenyon_common <- cbind(GO_Kenyon_common, pvalue_stall = GO_Kenyon_stall[i]$PValue)


GO_Hartl_common[, div := foldEnrich_stall / Fold.Enrichment]
GO_Kenyon_common[, div := foldEnrich_stall / Fold.Enrichment]

View(GO_Hartl_common)
View(GO_Kenyon_common)
write.csv(GO_Hartl_common, "/Users/KevinStein/Desktop/hartl.csv")
write.csv(GO_Kenyon_common, "/Users/KevinStein/Desktop/kenyon.csv")

GO_Hartl_subset <- GO_agg[Description == "tRNA aminoacylation for protein translation" |
                            Description == "germ cell development" |
                            Description == "proteasome complex" |
                            Description == "P granule" |
                            Description == "determination of adult lifespan" |
                            Description == "translation"]
GO_Hartl_subset <- GO_Hartl_subset[Dataset == "Hartl"]
i <- cbind(match(GO_Hartl_subset$Term, GO_Hartl_common$Term))
GO_Hartl_subset <- cbind(GO_Hartl_subset, div = GO_Hartl_common[i]$div)

test <- matrix(c(GO_Hartl_subset[Group == "stall" & Description == "tRNA aminoacylation for protein translation"]$Count,
                 GO_Hartl_subset[Group == "all" & Description == "tRNA aminoacylation for protein translation"]$Count,
                 GO_Hartl_subset[Group == "stall" & Description == "tRNA aminoacylation for protein translation"]$List.Total - 
                   GO_Hartl_subset[Group == "stall" & Description == "tRNA aminoacylation for protein translation"]$Count,
                 GO_Hartl_subset[Group == "all" & Description == "tRNA aminoacylation for protein translation"]$List.Total - 
                   GO_Hartl_subset[Group == "all" & Description == "tRNA aminoacylation for protein translation"]$Count), nrow = 2)
fisher.test(test)
# tRNA aminoacylation for protein translation: p-value = 0.06948, 2.33605
# germ cell development: p-value = 0.02372, odds = 2.125751
# proteasome complex: p-value = 0.1934, odds = 1.785842
# P granule: p-value = 0.1503, odds = 1.874646 
# determination of adult lifespan: p-value = 0.007698, odds = 1.743485
# translation: p-value = 0.2915, odds = 1.376072


GO_Kenyon_subset <- GO_agg[Description == "tRNA aminoacylation for protein translation" |
                             Description == "proteasome complex" |
                             Description == "translation" |
                             Description == "determination of adult lifespan" |
                             Description == "reproduction" |
                             Description == "cytosolic large ribosomal subunit"]
GO_Kenyon_subset <- GO_Kenyon_subset[Dataset == "Kenyon"]
i <- cbind(match(GO_Kenyon_subset$Term, GO_Kenyon_common$Term))
GO_Kenyon_subset <- cbind(GO_Kenyon_subset, div = GO_Kenyon_common[i]$div)

test <- matrix(c(GO_Kenyon_subset[Group == "stall" & Description == "translation"]$Count,
                 GO_Kenyon_subset[Group == "all" & Description == "translation"]$Count,
                 GO_Kenyon_subset[Group == "stall" & Description == "translation"]$List.Total - 
                   GO_Kenyon_subset[Group == "stall" & Description == "translation"]$Count,
                 GO_Kenyon_subset[Group == "all" & Description == "translation"]$List.Total - 
                   GO_Kenyon_subset[Group == "all" & Description == "translation"]$Count), nrow = 2)
fisher.test(test)
# tRNA aminoacylation for protein translation: p-value = 0.1834, odds = 1.869374
# proteasome complex: p-value = 0.5669, odds = 1.351985 
# translation: p-value = 0.3262, odds = 1.341277
# determination of adult lifespan: p-value = 0.2456, odds = 1.372844
# reproduction: p-value = 0.09665, odds = 1.518118
# cytosolic large ribosomal subunit: p-value = 0.3395, odds = 1.394318


plot <- ggplot(GO_Hartl_subset, aes(x = reorder(Description, div), y = Fold.Enrichment, fill = factor(Group))) + 
  geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("all","stall"), values = c("#80B1D3", "#FDB462"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = c(0.85,0.15), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_Hartl_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(GO_Kenyon_subset, aes(x = reorder(Description, div), y = Fold.Enrichment, fill = factor(Group))) + 
  geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("all","stall"), values = c("#FB9A99", "#E31A1C"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = c(0.85,0.15), axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_Kenyon_aggregates.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
