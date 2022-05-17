### Codon frequency in stalls ###
#set.seed(61)
#random_codon <- sample(sc_dtA[position > 20 & stopdist < 20]$codon, size = 1000000, replace = T)
proteome_codon <- as.data.table(table(sc_dtA[position > 20 & stopdist < -20 & residue != "X"]$codon))
codon_freq <- as.data.table(table(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$codon))
codon_freq <- cbind(codon_freq, proteome_codon[, 2])
setnames(codon_freq, c("codon", "stalls", "proteome"))
codon_freq[, stalls_freq := (stalls / sum(stalls))*100]
codon_freq[, proteome_freq := (proteome / sum(proteome))*100]
codon_freq[, div := stalls_freq / proteome_freq]
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(codon_freq$codon, codon_table$codon))
codon_freq <- cbind(codon_freq, residue = codon_table[i]$residue)

codon_freq$codon <- factor(codon_freq$codon, levels = c("GCA", "GCC", "GCG", "GCT",
                                                        "TGC", "TGT",
                                                        "GAC", "GAT",
                                                        "GAA", "GAG",
                                                        "TTC", "TTT",
                                                        "GGA", "GGC", "GGG", "GGT",
                                                        "CAC", "CAT",
                                                        "ATA", "ATC", "ATT",
                                                        "AAA", "AAG",
                                                        "CTA", "CTC", "CTG", "CTT", "TTA", "TTG",
                                                        "ATG",
                                                        "AAC", "AAT",
                                                        "CCA", "CCC", "CCG", "CCT",
                                                        "CAA", "CAG",
                                                        "CGA", "CGC", "CGG", "CGT", "AGA", "AGG",
                                                        "TCA", "TCC", "TCG", "TCT", "AGC", "AGT",
                                                        "ACA", "ACC", "ACG", "ACT",
                                                        "GTA", "GTC", "GTG", "GTT",
                                                        "TGG",
                                                        "TAA", "TAG", "TGA",
                                                        "TAC", "TAT"))

ggplot(codon_freq[codon != 'TAA' & codon != 'TAG' & codon != 'TGA'], aes(x = codon, y = log2(div), fill = residue)) + geom_col() +
  geom_hline(yintercept = 0, size = 0.5, color = 'black')
plot + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 5), axis.ticks = element_blank(), axis.line = element_blank())

i <- cbind(match(codon_freq$codon, codonusage_dt$codon))
codon_freq <- cbind(codon_freq, fraction = codonusage_dt[i]$fraction)
codon_freq <- cbind(codon_freq, per1000 = codonusage_dt[i]$per1000)
i <- cbind(match(codon_freq$codon, TE_dt$codon))
codon_freq <- cbind(codon_freq, norm.scale = TE_dt[i]$norm.scale)
codon_freq <- cbind(codon_freq, norm = TE_dt[i]$norm)
codon_freq <- cbind(codon_freq, classical.scale = TE_dt[i]$classical.scale)
codon_freq <- cbind(codon_freq, classical = TE_dt[i]$classical)
i <- cbind(match(codon_freq$codon, tAI_dt$codon))
codon_freq <- cbind(codon_freq, tAI = tAI_dt[i]$Scer)


ggplot(codon_freq, aes(fraction,div)) + geom_point() + ylab("pause site codon frequency") + xlab("codon usage (fraction)") 
ggplot(codon_freq, aes(classical.scale,div)) + geom_point() + ylab("pause site codon frequency") + xlab("codon usage (fraction)") 
ggplot(codon_freq, aes(tAI,div)) + geom_point() + ylab("pause site codon frequency") + xlab("tAI") 
ggplot(codon_freq[norm == "N"], aes("N",div)) + geom_boxplot() +
  geom_boxplot(data = codon_freq[norm == "O"], aes("O", div))
ggplot(codon_freq[norm == "N"], aes("N",div)) + geom_boxplot() +
  geom_boxplot(data = codon_freq[norm == "O"], aes("O", div)) + ylab("pause site codon frequency")


### Residue frequency in stalls ###
#set.seed(20)
#random_aa <- sample(sc_dtA[position > 20 & stopdist < -20]$residue, size = 100000, replace = T)
proteome_aa <- table(sc_dtA[position > 20 & stopdist < -20 & residue != "X"]$residue)
proteome_aa <- as.data.table(proteome_aa[c(1:19,21)])

WT_stalling_peaks_ageDep <- readRDS("WT_stalling_peaks_ageDep.rds")
aa_freq <- table(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$residue)
aa_freq <- as.data.table(aa_freq[c(1:19,21)])
aa_freq <- cbind(aa_freq, proteome_aa[, 2])
setnames(aa_freq, c("residue", "stalls", "proteome"))
aa_freq[, stalls_freq := (stalls / sum(stalls))*100]
aa_freq[, proteome_freq := (proteome / sum(proteome))*100]
aa_freq[, div := stalls_freq / proteome_freq]

WT_stalling_peaks_ageDep[, P := as.factor(substring(motif3,2,2))]
WT_stalling_peaks_ageDep[, E := as.factor(substring(motif3,1,1))]

aa_freqP <- data.table(table(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$P))
aa_freqP <- cbind(aa_freqP, proteome_aa[, 2])
setnames(aa_freqP, c("residue", "stalls", "proteome"))
aa_freqP[, stalls_freq := (stalls / sum(stalls))*100]
aa_freqP[, proteome_freq := (proteome / sum(proteome))*100]
aa_freqP[, div := stalls_freq / proteome_freq]

aa_freqE <- data.table(table(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$E))
aa_freqE <- cbind(aa_freqE, proteome_aa[, 2])
setnames(aa_freqE, c("residue", "stalls", "proteome"))
aa_freqE[, stalls_freq := (stalls / sum(stalls))*100]
aa_freqE[, proteome_freq := (proteome / sum(proteome))*100]
aa_freqE[, div := stalls_freq / proteome_freq]


ggplot(aa_freq[residue != 'X'], aes(residue, log2(div))) + geom_col()

test <- matrix(c(63,4988,(sum(aa_freq$stalls) - 63),(sum(aa_freq$proteome) - 4988)), nrow = 2)
fisher.test(test)



### Heat map ###
# x-axis is EPA
# y-axis is amino acid sorted by hydrophobicity
# fill is log2 enrichment


WT_stalling_peaks_ageDep <- readRDS("doc/WT_stalling_peaks_ageDep.rds")
sc_dtA <- readRDS("doc/sc_dtA.rds")
proteome_aa <- table(sc_dtA[position > 20 & stopdist < -20 & residue != "X"]$residue)
proteome_aa <- as.data.table(proteome_aa[c(1:19,21)])

aa_freq <- table(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$residue)
aa_freq <- as.data.table(aa_freq[c(1:19,21)])
aa_freq <- cbind(aa_freq, proteome_aa[, 2])
setnames(aa_freq, c("residue", "stalls", "proteome"))
aa_freq[, stalls_freq := (stalls / sum(stalls))*100]
aa_freq[, proteome_freq := (proteome / sum(proteome))*100]
aa_freq[, div := stalls_freq / proteome_freq]
aa_freq[, position := "3"]
aa_freq[, enrich := log2(div)]

WT_stalling_peaks_ageDep[, P := as.factor(substring(motif3,2,2))]
WT_stalling_peaks_ageDep[, E := as.factor(substring(motif3,1,1))]

aa_freqP <- data.table(table(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$P))
aa_freqP <- cbind(aa_freqP, proteome_aa[, 2])
setnames(aa_freqP, c("residue", "stalls", "proteome"))
aa_freqP[, stalls_freq := (stalls / sum(stalls))*100]
aa_freqP[, proteome_freq := (proteome / sum(proteome))*100]
aa_freqP[, div := stalls_freq / proteome_freq]
aa_freqP[, position := "2"]
aa_freqP[, enrich := log2(div)]

aa_freqE <- data.table(table(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$E))
aa_freqE <- cbind(aa_freqE, proteome_aa[, 2])
setnames(aa_freqE, c("residue", "stalls", "proteome"))
aa_freqE[, stalls_freq := (stalls / sum(stalls))*100]
aa_freqE[, proteome_freq := (proteome / sum(proteome))*100]
aa_freqE[, div := stalls_freq / proteome_freq]
aa_freqE[, position := "1"]
aa_freqE[, enrich := log2(div)]

aaHydro <- as.data.table(read.csv("/Users/KevinStein/Desktop/aaHydrophocity.csv"))

aa_freqfinal <- rbind(aa_freq, aa_freqP, aa_freqE)
i <- cbind(match(aa_freqfinal$residue, aaHydro$residue))
aa_freqfinal <- cbind(aa_freqfinal, hydro = aaHydro[i]$hydrophobicity)
aa_freqfinal <- cbind(aa_freqfinal, rank = aaHydro[i]$rank)


ggplot(aa_freqfinal, aes(position, rank, fill=enrich)) + geom_tile() + scale_fill_distiller(palette = "PuOr")
