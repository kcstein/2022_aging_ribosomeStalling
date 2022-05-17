### Codon frequency in stalls ###
#set.seed(61)
#random_codon <- sample(N2_dtA[position > 20 & stopdist < -20]$codon, size = 1000000, replace = F)
proteome_codon <- as.data.table(table(N2_dtA[position > 20 & stopdist < -20 & residue != "X"]$codon))
stalling_peaks_ageDep[, P := as.factor(substring(motif3,2,2))]
stalling_peaks_ageDep[, E := as.factor(substring(motif3,1,1))]

codon_freq <- as.data.table(table(stalling_peaks_ageDep[D12_pause > 10]$codon))
codon_freq <- cbind(codon_freq, proteome_codon[, 2])
setnames(codon_freq, c("codon", "stalls", "proteome"))
codon_freq[, stalls_freq := (stalls / sum(stalls))*100]
codon_freq[, proteome_freq := (proteome / sum(proteome))*100]
codon_freq[, div := stalls_freq / proteome_freq]
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(codon_freq$codon, codon_table$codon))
codon_freq <- cbind(codon_freq, residue = codon_table[i]$residue)

test <- matrix(c(18,11020,(sum(codon_freq$stalls) - 18),(sum(codon_freq$proteome) - 11020)), nrow = 2)
fisher.test(test)

i <- cbind(match(codon_freq$codon, codonusage_dt$codon))
codon_freq <- cbind(codon_freq, fraction = codonusage_dt[i]$fraction)
codon_freq <- cbind(codon_freq, per1000 = codonusage_dt[i]$per1000)

ggplot(codon_freq, aes(fraction,div)) + geom_point() + ylab("pause site codon frequency") + xlab("codon usage (fraction)") 
ggplot(codon_freq[fraction < 0.25], aes("N",div)) + geom_boxplot() +
  geom_boxplot(data = codon_freq[fraction > 0.25], aes("O", div))

wilcox.test(codon_freq[fraction < 0.25]$div, codon_freq[fraction > 0.25]$div)

# codon_freq <- data.table(stalls = table(stalling_peaks[D12_pause > 10]$codon))
# setnames(codon_freq, c("codon", "stalls"))
# i <- cbind(match(codon_freq$codon, codonusage_dt$codon))
# codon_freq <- cbind(codon_freq, proteome = codonusage_dt[i]$number)
# codon_freq[, stalls_freq := (stalls / sum(stalls))*100]
# codon_freq[, proteome_freq := (proteome / sum(proteome))*100]
# codon_freq[, div := stalls_freq / proteome_freq]
# i <- cbind(match(codon_freq$codon, codonusage_dt$codon))
# codon_freq <- cbind(codon_freq, residue = codonusage_dt[i]$residue)

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

plot <- ggplot(codon_freq[codon != 'TAA' & codon != 'TAG' & codon != 'TGA'], aes(codon, log2(div), fill = residue)) + geom_col()
plot + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 5), axis.ticks = element_blank(), axis.line = element_blank())


### Residue frequency in stalls ###
#set.seed(20)
#random_aa <- sample(N2_dtA[position > 20 & stopdist < -20]$residue, size = 100000, replace = T)
#proteome_aafreq <- c(A=6.325907,C=2.076066,D=5.304005,E=6.526124,F=4.787855,G=5.344524,H=2.281571,
 #       I=6.213785,K=6.339429,L=8.643033,M=2.648252,N=4.877476,P=4.892958,Q=4.081622,
  #      R=5.137781,S=8.080618,T=5.891986,V=6.229459,W=1.107328,Y=3.210222)
proteome_aa <- table(N2_dtA[position > 20 & stopdist < -20 & residue != "X"]$residue)
proteome_aa <- as.data.table(proteome_aa[c(1:19,21)])
stalling_peaks_ageDep[, P := as.factor(substring(motif3,2,2))]
stalling_peaks_ageDep[, E := as.factor(substring(motif3,1,1))]

stalling_peaks_ageDep <- readRDS("stalling_peaks_ageDep.rds")
aa_freq <- table(stalling_peaks_ageDep[D12_pause > 10]$residue)
aa_freq <- as.data.table(aa_freq[c(1:19,21)])
aa_freq <- cbind(aa_freq, proteome_aa[, 2])
setnames(aa_freq, c("residue", "stalls", "proteome"))
aa_freq[, stalls_freq := (stalls / sum(stalls))*100]
aa_freq[, proteome_freq := (proteome / sum(proteome))*100]
aa_freq[, div := stalls_freq / proteome_freq]

aa_freqP <- data.table(table(stalling_peaks_ageDep[D12_pause > 10]$P))
aa_freqP <- cbind(aa_freqP, proteome_aa[, 2])
setnames(aa_freqP, c("residue", "stalls", "proteome"))
aa_freqP[, stalls_freq := (stalls / sum(stalls))*100]
aa_freqP[, proteome_freq := (proteome / sum(proteome))*100]
aa_freqP[, div := stalls_freq / proteome_freq]

aa_freqE <- data.table(table(stalling_peaks_ageDep[D12_pause > 10]$E))
aa_freqE <- cbind(aa_freqE, proteome_aa[, 2])
setnames(aa_freqE, c("residue", "stalls", "proteome"))
aa_freqE[, stalls_freq := (stalls / sum(stalls))*100]
aa_freqE[, proteome_freq := (proteome / sum(proteome))*100]
aa_freqE[, div := stalls_freq / proteome_freq]


test <- matrix(c(aa_freq[residue == "R"]$stalls,
                 aa_freq[residue == "R"]$proteome,
                 (sum(aa_freq$stalls) - aa_freq[residue == "R"]$stalls),
                 (sum(aa_freq$proteome) - aa_freq[residue == "R"]$proteome)), nrow = 2)
fisher.test(test)

ggplot(aa_freq[residue != 'X'], aes(residue, log2(div))) + geom_col()