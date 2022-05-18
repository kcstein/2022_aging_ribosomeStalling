### Parse data from Hartl and Kenyon papers to get proteins with age-dependent aggregation ###

library(data.table)
library(plyr)
library(dplyr)
library(tidyr)


agg.propensity <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Datasets/aggregationPropensity_WaltherHartl_Cell2015_TableS1E.csv", header = TRUE, stringsAsFactors = TRUE))
aggregates_hartl <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Datasets/D1_D12aggregation_WaltherHartl_Cell2015_TableS1D.csv", header = TRUE, stringsAsFactors = TRUE))
aggregates_kenyon <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Datasets/insolubleProteins_DavidKenyon_PLOSBio2010_TableS1.csv", header = TRUE, stringsAsFactors = TRUE))
proteome_hartl <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Datasets/totalProteome_WaltherHartl_Cell2015_TableS1B.csv", header = TRUE, stringsAsFactors = TRUE))

### Find proteins in Hartl dataset showing age-dependent aggregation ###

# "We isolated insoluble proteins from total lysates of WT animals by centrifugation 
# and performed MS analysis using lysate from labeled worms for quantification (Figure S4C). 
# About 90% of the proteins that were quantified in three out of four experiments 
# (975 of 1,083 proteins) accumulated significantly in the insoluble fraction of 
# day 12 animals relative to day 1 (Table S1D)."

for (i in 1:nrow(aggregates_hartl)) {
  D1_silac <- c(aggregates_hartl[i, c(3:6)])
  D12_silac <- c(aggregates_hartl[i, c(7:10)])
  aggregates_hartl[i, D1_silac_samples := length(na.omit(unlist(D1_silac)))]
  aggregates_hartl[i, D12_silac_samples := length(na.omit(unlist(D12_silac)))]
  aggregates_hartl[i, D1_silac_median := median(na.omit(unlist(D1_silac)))]
  aggregates_hartl[i, D12_silac_median := median(na.omit(unlist(D12_silac)))]
}
aggregates_hartl[, silac_median_diff := D12_silac_median - D1_silac_median]

# isolate proteins that have silac_diff > 0 and statsistically significant, 3 replicates
aggregates_hartl_silac <- aggregates_hartl[D1_silac_samples >= 3 & D12_silac_samples >= 3]
for (i in 1:nrow(aggregates_hartl_silac)) {
  D1_silac <- c(aggregates_hartl_silac[i, c(3:6)])
  D12_silac <- c(aggregates_hartl_silac[i, c(7:10)])
  silac_test <- t.test((na.omit(unlist(D1_silac))), (na.omit(unlist(D12_silac))))
  aggregates_hartl_silac[i, silac_pvalue := silac_test$p.value]
}
aggregates_hartl_silac[, silac_padj := p.adjust(silac_pvalue, method = "BH")]

length(aggregates_hartl_silac[silac_pvalue < 0.05 & silac_median_diff > 0]$Log2.Ratio.L.H.AgWTD01.1)
length(aggregates_hartl_silac[silac_padj < 0.05 & silac_median_diff > 0]$Log2.Ratio.L.H.AgWTD01.1)
aggregates_hartl_silac <- separate(aggregates_hartl_silac, col = Uniprot.ID, into = c("Uniprot.ID1", NA), sep = ";", remove = FALSE)


### Isolate proteins that change abundance with age ### 
# only those proteins were displayed that were quantified at day 1 and 
# at least at three consecutive time points, but since I'm only using data from Day 1 and 12, I will count as proteome the proteins that were quantified at these two points
proteome_hartl_subset <- proteome_hartl[!is.na(Log2.Ratio.L.H.NormalizedWTday01) & !is.na(Log2.Ratio.L.H.NormalizedWTday12)]
proteome_hartl_subset[, D1D12_diff := Log2.Ratio.L.H.NormalizedWTday12 - Log2.Ratio.L.H.NormalizedWTday01]
length(proteome_hartl_subset[D1D12_diff < 0]$Uniprot.ID)
proteome_hartl_subset <- separate(proteome_hartl_subset, col = Uniprot.ID, into = c("Uniprot.ID1", NA), sep = ";", remove = FALSE)
agg.propensity <- separate(agg.propensity, col = Uniprot.ID, into = c("Uniprot.ID1", NA), sep = ";", remove = FALSE)


### For downstream analysis: obtain sequences then analyze aa composition, GO, number proteins with KR stretches
write.table(aggregates_kenyon[Reproducible_hits == 1, 1], "aggregates_kenyon_uniprotID.csv", sep = ',', row.names = FALSE, col.names = FALSE)
write.table(aggregates_hartl_silac[silac_padj < 0.05 & silac_median_diff > 0, 2],
            "silac_sig_padj.csv", sep = ',', row.names = FALSE, col.names = FALSE)
write.table(agg.propensity, "hartl_proteome.csv", sep = ',', row.names = FALSE, col.names = FALSE)
write.table(proteome_hartl_subset[D1D12_diff > 0, 2], 
            "proteome_up.csv", sep = ',', row.names = FALSE, col.names = FALSE)
write.table(proteome_hartl_subset[D1D12_diff < 0, 2], 
            "proteome_down.csv", sep = ',', row.names = FALSE, col.names = FALSE)

# after obtaining gene names from uniprot, combine into one datatable for reference
hartl_genename <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_genename.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(hartl_genename) <- c("uniprot", "gene")
hartl_genename_unmapped <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_genename_unmapped.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(hartl_genename_unmapped) <- c("uniprot")
hartl_WBgene <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_WBgene.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(hartl_WBgene) <- c("uniprot", "WBgene")
hartl_WBgene_unmapped <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_WBgene_unmapped.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(hartl_WBgene_unmapped) <- c("uniprot")
hartl_WBtranscript <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_WBtranscript.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(hartl_WBtranscript) <- c("uniprot", "transcript")
hartl_WBtranscript_unmapped <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_WBtranscript_unmapped.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(hartl_WBtranscript_unmapped) <- c("uniprot")
kenyonAggregates_genename <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/kenyonAggregates_genename.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(kenyonAggregates_genename) <- c("uniprot", "gene")
kenyonAggregates_genename_unmapped <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/kenyonAggregates_genename_unmapped.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(kenyonAggregates_genename_unmapped) <- c("uniprot")
kenyonAggregates_WBgene <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/kenyonAggregates_WBgene.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(kenyonAggregates_WBgene) <- c("uniprot", "WBgene")
kenyonAggregates_WBgene_unmapped <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/kenyonAggregates_WBgene_unmapped.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(kenyonAggregates_WBgene_unmapped) <- c("uniprot")
kenyonAggregates_WBtranscript <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/kenyonAggregates_WBtranscript.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(kenyonAggregates_WBtranscript) <- c("uniprot", "transcript")
kenyonAggregates_WBtranscript_unmapped <- as.data.table(read.delim("/Users/KevinStein/Desktop/AggregatesAnalysis/kenyonAggregates_WBtranscript_unmapped.txt", header = TRUE, stringsAsFactors = TRUE))
colnames(kenyonAggregates_WBtranscript_unmapped) <- c("uniprot")

uniprot <- unique(hartl_WBgene$uniprot)
hartl_proteome <- data.table(uniprot = uniprot)
i <- cbind(match(hartl_proteome$uniprot, hartl_WBgene$uniprot))
hartl_proteome <- cbind(hartl_proteome, WBgene = hartl_WBgene[i]$WBgene)
i <- cbind(match(hartl_proteome$uniprot, hartl_WBtranscript$uniprot))
hartl_proteome <- cbind(hartl_proteome, WBtranscript = hartl_WBtranscript[i]$transcript)
i <- cbind(match(hartl_proteome$uniprot, hartl_genename$uniprot))
hartl_proteome <- cbind(hartl_proteome, gene = hartl_genename[i]$gene)
uniprot <- unique(hartl_WBgene_unmapped$uniprot)
hartl_proteome_unmapped <- data.table(uniprot = uniprot)
i <- cbind(match(hartl_proteome_unmapped$uniprot, hartl_WBgene_unmapped$uniprot))
hartl_proteome_unmapped <- cbind(hartl_proteome_unmapped, WBgene = hartl_WBgene_unmapped[i]$WBgene)
i <- cbind(match(hartl_proteome_unmapped$uniprot, hartl_WBtranscript_unmapped$uniprot))
hartl_proteome_unmapped <- cbind(hartl_proteome_unmapped, WBtranscript = hartl_WBtranscript_unmapped[i]$transcript)
i <- cbind(match(hartl_proteome_unmapped$uniprot, hartl_genename_unmapped$uniprot))
hartl_proteome_unmapped <- cbind(hartl_proteome_unmapped, gene = hartl_genename_unmapped[i]$gene)
uniprot <- unique(kenyonAggregates_WBgene$uniprot)
kenyonAggregates <- data.table(uniprot = uniprot)
i <- cbind(match(kenyonAggregates$uniprot, kenyonAggregates_WBgene$uniprot))
kenyonAggregates <- cbind(kenyonAggregates, WBgene = kenyonAggregates_WBgene[i]$WBgene)
i <- cbind(match(kenyonAggregates$uniprot, kenyonAggregates_WBtranscript$uniprot))
kenyonAggregates <- cbind(kenyonAggregates, WBtranscript = kenyonAggregates_WBtranscript[i]$transcript)
i <- cbind(match(kenyonAggregates$uniprot, kenyonAggregates_genename$uniprot))
kenyonAggregates <- cbind(kenyonAggregates, gene = kenyonAggregates_genename[i]$gene)
uniprot <- unique(kenyonAggregates_WBgene_unmapped$uniprot)
kenyonAggregates_unmapped <- data.table(uniprot = uniprot)
i <- cbind(match(kenyonAggregates_unmapped$uniprot, kenyonAggregates_WBgene_unmapped$uniprot))
kenyonAggregates_unmapped <- cbind(kenyonAggregates_unmapped, WBgene = kenyonAggregates_WBgene_unmapped[i]$WBgene)
i <- cbind(match(kenyonAggregates_unmapped$uniprot, kenyonAggregates_WBtranscript_unmapped$uniprot))
kenyonAggregates_unmapped <- cbind(kenyonAggregates_unmapped, WBtranscript = kenyonAggregates_WBtranscript_unmapped[i]$transcript)
i <- cbind(match(kenyonAggregates_unmapped$uniprot, kenyonAggregates_genename_unmapped$uniprot))
kenyonAggregates_unmapped <- cbind(kenyonAggregates_unmapped, gene = kenyonAggregates_genename_unmapped[i]$gene)

hartl_aggregates <- as.data.table(read.csv("/Users/KevinStein/Desktop/AggregatesAnalysis/proteins/silac_sig_padj.csv", header = FALSE, stringsAsFactors = TRUE))
hartl_aggregates <- hartl_proteome[hartl_proteome$uniprot %in% hartl_aggregates$V1]
hartl_proteome_up <- as.data.table(read.csv("/Users/KevinStein/Desktop/AggregatesAnalysis/proteins/proteome_up.csv", header = FALSE, stringsAsFactors = TRUE))
hartl_proteome_up <- hartl_proteome[hartl_proteome$uniprot %in% hartl_proteome_up$V1]
hartl_proteome_down <- as.data.table(read.csv("/Users/KevinStein/Desktop/AggregatesAnalysis/proteins/proteome_down.csv", header = FALSE, stringsAsFactors = TRUE))
hartl_proteome_down <- hartl_proteome[hartl_proteome$uniprot %in% hartl_proteome_down$V1]

write.csv(hartl_aggregates, "/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_aggregates.csv", row.names = FALSE)
write.csv(hartl_proteome_up, "/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_proteome_up.csv", row.names = FALSE)
write.csv(hartl_proteome_down, "/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_proteome_down.csv", row.names = FALSE)
write.csv(hartl_proteome, "/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_proteome.csv", row.names = FALSE)
write.csv(hartl_proteome_unmapped, "/Users/KevinStein/Desktop/AggregatesAnalysis/hartl_proteome_unmapped.csv", row.names = FALSE)
write.csv(kenyonAggregates, "/Users/KevinStein/Desktop/AggregatesAnalysis/kenyon_aggregates.csv", row.names = FALSE)
write.csv(kenyonAggregates_unmapped, "/Users/KevinStein/Desktop/AggregatesAnalysis/kenyon_aggregates_unmapped.csv", row.names = FALSE)

# Add gene name info to agg.propensity dataset
i <- cbind(match(agg.propensity$Uniprot.ID1, hartl_proteome$uniprot))
agg.propensity <- cbind(agg.propensity, WBgene = hartl_proteome[i]$WBgene)
agg.propensity <- cbind(agg.propensity, WBtranscript = hartl_proteome[i]$WBtranscript)
agg.propensity <- cbind(agg.propensity, gene = hartl_proteome[i]$gene)
write.csv(agg.propensity, "/Users/KevinStein/Desktop/AggregatesAnalysis/agg.propensity.csv", row.names = FALSE)