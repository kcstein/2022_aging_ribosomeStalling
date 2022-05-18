library(data.table)
library(plyr)
library(dplyr)
library(tidyr)


### Correlate stalling and aggregation ###
stalling_peaks_ageDep <- readRDS("stalling_peaks_ageDep.rds")
hartl_aggregates <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/proteinLists/hartl_aggregates.csv', header = FALSE, stringsAsFactors = TRUE))
colnames(hartl_aggregates) <- c("orf", "WBgene")
kenyon_aggregates <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/proteinLists/kenyon_aggregates.csv', header = FALSE, stringsAsFactors = TRUE))
colnames(kenyon_aggregates) <- c("orf", "WBgene")
hartl_proteome <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/proteinLists/hartl_proteome.csv', header = FALSE, stringsAsFactors = TRUE))
colnames(hartl_proteome) <- c("orf", "WBgene")
genenames <- as.data.table(read.delim('/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/WBcel235.Ensembl95/Cele_proteincoding_cdsLength.txt', header = TRUE, stringsAsFactors = TRUE))
i <- cbind(match(stalling_peaks_ageDep$orf, genenames$Transcript.stable.ID))
stalling_peaks_ageDep <- cbind(stalling_peaks_ageDep, WBgene = genenames[i]$Gene.stable.ID)
i <- cbind(match(stalling_peaks_ageDep$WBgene, hartl_proteome$WBgene))
stalling_peaks_ageDep <- cbind(stalling_peaks_ageDep, gene = hartl_proteome[i]$gene)
stalling_peaks_ageDep <- cbind(stalling_peaks_ageDep, uniprot = hartl_proteome[i]$uniprot)
stalling_peaks_orfs <- stalling_peaks_ageDep[, .SD[which.max(peak)], by = orf]
#saveRDS(stalling_peaks_ageDep, "stalling_peaks_ageDep.rds")

N2_dtA <- readRDS("N2_dtA.rds")
orfs <- N2_dtA[, .SD[which.max(position)], by = orf]
i <- cbind(match(orfs$orf, genenames$Transcript.stable.ID))
orfs <- cbind(orfs, WBgene = genenames[i]$Gene.stable.ID)
i <- cbind(match(N2_dtA$orf, genenames$Transcript.stable.ID))
N2_dtA <- cbind(N2_dtA, WBgene = genenames[i]$Gene.stable.ID)
i <- cbind(match(N2_dtA$WBgene, hartl_proteome$WBgene))
N2_dtA <- cbind(N2_dtA, gene = hartl_proteome[i]$gene)
N2_dtA <- cbind(N2_dtA, uniprot = hartl_proteome[i]$uniprot)
#saveRDS(N2_dtA, "N2_dtA.rds")


# proteome: 6961 (no aggregate and no stall: 5651 (6961 - 742 - 730 + 162))
hartl_proteome_orfs <- as.data.table(unique(hartl_proteome$WBgene))
length(hartl_proteome_orfs$V1)
setnames(hartl_proteome_orfs, c("WBgene"))
i <- cbind(match(hartl_proteome_orfs$WBgene, orfs$WBgene))
hartl_proteome_orfs <- cbind(hartl_proteome_orfs, orf = orfs[i]$orf)
write.csv(hartl_proteome_orfs, "hartl_proteome.csv")


# stall: 742 (stall no aggregate = 742-162)
length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene)
write.csv(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$V1], "stalling_peaks_hartlProteome.csv")

# aggregate: 730 (aggregate no stall = 734-162)
hartl_aggregates_orfs <- as.data.table(unique(hartl_aggregates$WBgene))
length(hartl_aggregates_orfs[hartl_aggregates_orfs$V1 %in% hartl_proteome_orfs$WBgene]$V1)
setnames(hartl_aggregates_orfs, c("WBgene"))
i <- cbind(match(hartl_aggregates_orfs$WBgene, orfs$WBgene))
hartl_aggregates_orfs <- cbind(hartl_aggregates_orfs, orf = orfs[i]$orf)
write.csv(hartl_aggregates_orfs, "hartl_aggregates.csv")

# aggregate and stall: 162
length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)
write.csv(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene], "hartl_aggregatesANDstall.csv")
write.csv(hartl_aggregates_orfs[!hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene], "hartl_aggregatesNOTstall.csv")
write.csv(stalling_peaks_orfs[!stalling_peaks_orfs$WBgene %in% hartl_aggregates_orfs$WBgene], "stallNOThartl_aggregates.csv")

# kenyon aggregate: 384 (aggregate no stall = 384-109)
kenyon_aggregates_orfs <- as.data.table(unique(kenyon_aggregates$WBgene))
length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$V1 %in% hartl_proteome_orfs$WBgene]$V1)
setnames(kenyon_aggregates_orfs, c("WBgene"))
i <- cbind(match(kenyon_aggregates_orfs$WBgene, orfs$WBgene))
kenyon_aggregates_orfs <- cbind(kenyon_aggregates_orfs, orf = orfs[i]$orf)
write.csv(kenyon_aggregates_orfs, "kenyon_aggregates.csv")

# kenyon aggregate and stall: 109
length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)
write.csv(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene], "kenyon_aggregatesANDstall.csv")
write.csv(kenyon_aggregates_orfs[!kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene], "kenyon_aggregatesNOTstall.csv")
write.csv(stalling_peaks_orfs[!stalling_peaks_orfs$WBgene %in% kenyon_aggregates_orfs$WBgene], "stallNOTkenyon_aggregates.csv")

test <- matrix(c(length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene), # aggregation + stall
                 length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                   length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene),
                 length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                   length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene),
                 length(hartl_proteome_orfs$WBgene) - 
                   length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                   length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) +
                   length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)), nrow = 2)
test1 <- fisher.test(test)
test1$p.value

test <- matrix(c(length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene),
                 length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                   length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene),
                 length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                   length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene),
                 length(hartl_proteome_orfs$WBgene) - 
                   length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) -
                   length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) +
                   length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)), nrow = 2)
test1 <- fisher.test(test)
test1$p.value


# bar chart: fraction of datasets that aggregate, i.e. fraction of proteins that aggregate in 
# the proteome (730 / 6961) VS fraction of proteins that aggregate in the stalling pool (162 / 742)

temp <- data.table(Category = c("Proteome", "Stalling"), 
                   FractionAggregation = c((length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) /
                                              length(hartl_proteome_orfs$WBgene)),
                                           (length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)) /
                                             length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene)))

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

temp <- data.table(Category = c("Proteome", "Stalling"), 
                   FractionAggregation = c((length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene) /
                                              length(hartl_proteome_orfs$WBgene)),
                                           (length(hartl_aggregates_orfs[hartl_aggregates_orfs$WBgene %in% stalling_peaks_orfs$WBgene]$WBgene)) /
                                             length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$WBgene]$WBgene)))
plot <- ggplot(temp, aes(x = Category, y=FractionAggregation, fill = Category)) + geom_col(color = "black", size = 0.5) + scale_y_continuous(limits = c(0,0.25), expand = c(0.0008,0.0008)) +
  scale_fill_manual(limits = c("Proteome","Stalling"), values = c("#80B1D3", "#FDB462"), name = "")
plot <- plot + theme_classic(18) + labs(y = "Fraction of aggregated\nproteins in dataset", x = "") +
  theme(legend.position = "none", legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Stalling_HartlAggregation_p5.6e-22.pdf", plot, width = 3.2, height = 4, dpi = 300, useDingbats = F)


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


temp <- data.table(Category = c("Proteome", "Stalling"), 
                   FractionAggregation = c((length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$V1 %in% hartl_proteome_orfs$V1]$V1) /
                                              length(hartl_proteome_orfs$V1)),
                                           (length(kenyon_aggregates_orfs[kenyon_aggregates_orfs$V1 %in% stalling_peaks_orfs$WBgene]$V1)) /
                                             length(stalling_peaks_orfs[stalling_peaks_orfs$WBgene %in% hartl_proteome_orfs$V1]$WBgene)))
plot <- ggplot(temp, aes(x = Category, y=FractionAggregation, fill = Category)) + geom_col(color = "black", size = 0.5) + scale_y_continuous(limits = c(0,0.15), expand = c(0.0005,0.0005)) +
  scale_fill_manual(limits = c("Proteome","Stalling"), values = c("gray50", "#CB181D"), name = "")
plot <- plot + theme_classic(18) + labs(y = "Fraction of aggregated\nproteins in dataset", x = "") +
  theme(legend.position = "none", legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Stalling_KenyonAggregation_p1.6e-23.pdf", plot, width = 3.2, height = 4, dpi = 300, useDingbats = F)
