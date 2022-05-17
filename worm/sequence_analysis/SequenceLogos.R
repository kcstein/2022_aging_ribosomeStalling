### References
# http://bioconductor.org/packages/devel/bioc/vignettes/Logolas/inst/doc/Logolas.html
# https://omarwagih.github.io/ggseqlogo/
# http://kplogo.wi.mit.edu/index.html?
# https://dongyuexie.github.io/Sequence-Motif-Depletion/index.html
# https://iomics.ugent.be/icelogoserver/

library(data.table)
library(ggplot2)
#devtools::install_github("collectivemedia/tictoc")
library(tictoc)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')

### Logo packages ###
install.packages("devtools")
library(devtools)
install.packages("RWebLogo")
library(RWebLogo)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Logolas")
library(Logolas)


### Subset sequences ###
stalling_peaks_seq <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/Worm/doc/Fishers/stalling_peaks_seq.csv", header = TRUE, stringsAsFactors = FALSE))
stalling_peaks_seq60 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/Worm/doc/Fishers/stalling_peaks_60mer_seq.csv", header = TRUE, stringsAsFactors = FALSE))
stalling_peaks_seq[, tunnel := substring(sequence,1,27)]
stalling_peaks_seq[, active := substring(sequence,28,30)]
stalling_peaks_seq60[, nascent := substring(sequence,1,29)]
stalling_peaks_seq60[, tunnel := substring(sequence,30,57)]
stalling_peaks_seq60[, active := substring(sequence,58,60)]


### Logos ###
# for use with kplogo: A:0.06325907,C:0.02076066,D:0.05304005,E:0.06526124,F:0.04787855,G:0.05344524,H:0.02281571,I:0.06213785,K:0.06339429,L:0.08643033,M:0.02648252,N:0.04877476,P:0.04892958,Q:0.04081622,R:0.05137781,S:0.08080618,T:0.05891986,V:0.06229459,W:0.01107328,Y:0.03210222
stalling_peaks_seq <- readRDS("stalling_peaks_seq.rds")
#bg <- c(A=0.06325907,C=0.02076066,D=0.05304005,E=0.06526124,F=0.04787855,G=0.05344524,H=0.02281571,I=0.06213785,K=0.06339429,L=0.08643033,M=0.02648252,N=0.04877476,P=0.04892958,Q=0.04081622,R=0.05137781,S=0.08080618,T=0.05891986,V=0.06229459,W=0.01107328,Y=0.03210222)
bg <- c(A=0.06337123,C=0.02081661,D=0.05376245,E=0.06618098,F=0.04724033,G=0.05420618,
        H=0.02301501,I=0.06203379,K=0.06325141,L=0.08587701,M=0.02462631,N=0.04885255,
        P=0.04920322,Q=0.04115967,R=0.05131655,S=0.07982485,T=0.05913152,V=0.06248138,W=0.01124630,Y=0.03240266)

weblogo(stalling_peaks_seq[D12_pause > 10]$tunnel, open = TRUE, file.out = "stalling_fishers_pause10.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
weblogo(stalling_peaks_seq[D12_pause > 10]$tunnel, open = TRUE, file.out = "stalling_fishers_pause10IC.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        composition = c('A'=6.325907,'C'=2.076066,'D'=5.304005,'E'=6.526124,'F'=4.787855,'G'=5.344524,'H'=2.281571,'I'=6.213785,'K'=6.339429,'L'=8.643033,'M'=2.648252,'N'=4.877476,'P'=4.892958,'Q'=4.081622,'R'=5.137781,'S'=8.080618,'T'=5.891986,'V'=6.229459,'W'=1.107328,'Y'=3.210222),
        color.scheme = 'chemistry', errorbars = FALSE, yaxis = 0.05)
logomaker(stalling_peaks_seq$active, type = "EDLogo", bg = bg)
logomaker(stalling_peaks_seq60[D12_pause > 6]$nascent, type = "EDLogo", bg = bg, color_seed = 6)

N2_fishers <- readRDS("N2_fishers.rds")
bg <- summary(N2_fishers[position > 20 & stopdist < -20]$residue) / 
  length(N2_fishers[position > 20 & stopdist < -20]$residue)
bg <- bg[c(1:19,21)]

logomaker(stalling_peaks_ageDep[D12_pause > 10 & residue != "W"]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
summary(stalling_peaks_ageDep$residue) # Trp under 100

# weblogo(stalling_peaks_ageDep[D12_pause > 10]$motif3, open = TRUE, file.out = "stalling_fishers_pause10.pdf",
#         format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
#         color.scheme = 'chemistry', units = 'probability')
# weblogo(stalling_peaks_ageDep[D12_pause > 10]$motif3, open = TRUE, file.out = "stalling_fishers_pause10IC.pdf",
#         format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
#         composition = c('A'=6.325907,'C'=2.076066,'D'=5.304005,'E'=6.526124,'F'=4.787855,'G'=5.344524,'H'=2.281571,'I'=6.213785,'K'=6.339429,'L'=8.643033,'M'=2.648252,'N'=4.877476,'P'=4.892958,'Q'=4.081622,'R'=5.137781,'S'=8.080618,'T'=5.891986,'V'=6.229459,'W'=1.107328,'Y'=3.210222),
#         color.scheme = 'chemistry', errorbars = FALSE, yaxis = 0.05)
summary(stalling_peaks_ageDep[D12_pause > 10]$residue) / 
  length(stalling_peaks_ageDep[D12_pause > 10]$residue)


### Genome-wide tripeptide active site analysis ###
N2_dtA <- readRDS("N2_dtA.rds")
temp <- N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & 
                 N2_D6_2A_rpc >= 0.5 & 
                 D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 &
                 D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                 D12_1A_sum >= 64 & D12_2A_sum >= 64 &
                 N2_D6_2A_sum >= 64 & 
                 D1_1A > 0 & !is.na(D1_1A_pause) & D1_2A > 0 & !is.na(D1_2A_pause) &
                 N2_D6_2A > 0 & !is.na(N2_D6_2A_pause) &
                 D12_1A > 0 & !is.na(D12_1A_pause) & D12_2A > 0 & !is.na(D12_2A_pause) &
                 position > 20 & stopdist < -20]
motifs <- unique(temp[!is.na(motif3)]$motif3)

tic("runtime")
motif <- NULL
count <- NULL
D1 <- NULL
D1_med <- NULL
D6 <- NULL
D6_med <- NULL
D12 <- NULL
D12_med <- NULL
p_12_1 <- NULL
p_6_1 <- NULL
p_12_6 <- NULL
for (i in motifs) {
  print(match(i, motifs))
  motif <- c(motif, i)
  count1 <- length(which(temp$motif3 == i))
  count <- c(count, count1)
  D1.1 <- mean(temp[motif3 == i]$D1_pause)
  D1 <- c(D1, D1.1)
  D1.2 <- median(temp[motif3 == i]$D1_pause)
  D1_med <- c(D1_med, D1.2)
  D6.1 <- mean(temp[motif3 == i]$N2_D6_2A_pause)
  D6 <- c(D6, D6.1)
  D6.2 <- median(temp[motif3 == i]$N2_D6_2A_pause)
  D6_med <- c(D6_med, D6.2)
  D12.1 <- mean(temp[motif3 == i]$D12_pause)
  D12 <- c(D12, D12.1)
  D12.2 <- median(temp[motif3 == i]$D12_pause)
  D12_med <- c(D12_med, D12.2)
  test_12_1 <- wilcox.test(temp[motif3 == i]$D1_pause, temp[motif3 == i]$D12_pause)
  p12.1 <- test_12_1$p.value
  p_12_1 <- c(p_12_1, p12.1)
  test_6_1 <- wilcox.test(temp[motif3 == i]$D1_pause, temp[motif3 == i]$N2_D6_2A_pause)
  p6.1 <- test_6_1$p.value
  p_6_1 <- c(p_6_1, p6.1)
  test_12_6 <- wilcox.test(temp[motif3 == i]$N2_D6_2A_pause, temp[motif3 == i]$D12_pause)
  p12.6 <- test_12_6$p.value
  p_12_6 <- c(p_12_6, p12.6)
}
toc()
motif_dt <- data.table(motif = motif, count = count, 
                       D1_pause = D1, D6_pause = D6, D12_pause = D12,
                       D1_pause_med = D1_med, D6_pause_med = D6_med, D12_pause_med = D12_med, 
                       p_12_1 = p_12_1, p_6_1 = p_6_1, p_12_6 = p_12_6)
motif_dt[, ratio_12_1 := D12_pause / D1_pause]
motif_dt[, ratio_med_12_1 := D12_pause_med / D1_pause_med]
motif_dt[, padj_12_1 := p.adjust(p_12_1, method = "BH")]
motif_dt[, ratio_6_1 := D6_pause / D1_pause]
motif_dt[, ratio_med_6_1 := D6_pause_med / D1_pause_med]
motif_dt[, padj_6_1 := p.adjust(p_6_1, method = "BH")]
motif_dt[, ratio_12_6 := D12_pause / D6_pause]
motif_dt[, ratio_med_12_6 := D12_pause_med / D6_pause_med]
motif_dt[, padj_12_6 := p.adjust(p_12_6, method = "BH")]
#saveRDS(motif_dt, "motif_dt.rds")

ggplot(motif_dt[count >= 100], aes(D1_pause, D12_pause)) + geom_point() + xlim(1,2.5) + ylim(1,2.5) + 
  geom_point(data = motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2], aes(D1_pause, D12_pause), color = "blue") + 
  geom_point(data = motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2 & p_12_1 < 0.05], aes(D1_pause, D12_pause), color = "red") + 
  geom_abline(slope = 1, intercept = 0, color = 'red') 

# count >= 100 is important; median and mean don't matter; all examples with p < 0.05 also have D12_pause > 1
# for top 50 or 100, median is better
# selecting for highest D12_pause is biasing to proline

motif_dt <- readRDS("motif_dt.rds")
weblogo(motif_dt[count >= 100 & ratio > 1 & D12_pause > 2]$motif, open = TRUE, file.out = "winner_Count_D12above1.5_ratio1.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
weblogo(motif_dt[count >= 100 & ratio > 1 & D12_pause > 2]$motif, open = TRUE, file.out = "Count_p_IC.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        composition = c('A'=6.325907,'C'=2.076066,'D'=5.304005,'E'=6.526124,'F'=4.787855,'G'=5.344524,'H'=2.281571,'I'=6.213785,'K'=6.339429,'L'=8.643033,'M'=2.648252,'N'=4.877476,'P'=4.892958,'Q'=4.081622,'R'=5.137781,'S'=8.080618,'T'=5.891986,'V'=6.229459,'W'=1.107328,'Y'=3.210222),
        color.scheme = 'chemistry', errorbars = FALSE, yaxis = 1)
# top 100 ranked by D12_pause
weblogo(test3$motif, open = TRUE, file.out = "test1.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
#logomaker(motif_dt[count >= 100 & ratio > 1 & pvalue < 0.05 & D12_pause > 1.5]$motif, type = "EDLogo", bg = bg)
logomaker(motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_6_1 > 1 & D6_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2]$motif, type = "Logo", bg = bg, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_6_1 > 1 & D6_pause > 2]$motif, type = "Logo", bg = bg, color_seed = 6)

bg_noW <- bg[c(1:18,20)]
logomaker(motif_dt[count >= 100 & ratio_12_6 > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg_noW, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio_12_6 > 1 & D12_pause > 2]$motif, type = "Logo", bg = bg_noW, color_seed = 6)

motif_dt[, E := as.factor(substring(motif,1,1))]
motif_dt[, P := as.factor(substring(motif,2,2))]
motif_dt[, A := as.factor(substring(motif,3,3))]
summary(motif_dt[count >= 100 & ratio_12_6 > 1]$E)
summary(motif_dt[count >= 100 & ratio_12_6 > 1]$P)
summary(motif_dt[count >= 100 & ratio_12_6 > 1]$A)

logomaker(motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2 & p_12_1 < 0.05]$motif, type = "EDLogo", bg = bg, color_seed = 6)
motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2 & p_12_1 < 0.05]$motif

ggplot(motif_dt[count >= 100], aes(log2(ratio_12_1), D12_pause)) + geom_point() +
  geom_point(data = motif_dt[count >= 100 & E == "R"], aes(log2(ratio_12_1), D12_pause), color = "blue") +
  geom_point(data = motif_dt[count >= 100 & P == "R"], aes(log2(ratio_12_1), D12_pause), color = "orange") +
  geom_point(data = motif_dt[count >= 100 & A == "R"], aes(log2(ratio_12_1), D12_pause), color = "red")
  
ggplot(motif_dt[count >= 100], aes(log2(ratio_12_1), D12_pause)) + geom_point() +
  geom_point(data = motif_dt[count >= 100 & ratio_12_1 > 1 & D12_pause > 2 & p_12_1 < 0.05], aes(log2(ratio_12_1), D12_pause), color = 'red')

# to use with kplogo online server
write.csv(motif_dt[count >= 100 & ratio > 1 & D12_pause > 2], "/Users/KevinStein/Desktop/worm_ratio1_D12_1.5.csv")


### Genome-wide dipeptide motif analysis ###
temp <- N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & 
                 D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 &
                 D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                 D12_1A_sum >= 64 & D12_2A_sum >= 64 &
                 D1_1A > 0 & !is.na(D1_1A_pause) & D1_2A > 0 & !is.na(D1_2A_pause) &
                 D12_1A > 0 & !is.na(D12_1A_pause) & D12_2A > 0 & !is.na(D12_2A_pause) &
                 position > 15 & stopdist < -15]
motifs <- unique(temp[!is.na(motif2)]$motif2)

tic("runtime")
motif <- NULL
count <- NULL
D1 <- NULL
D1_med <- NULL
D12 <- NULL
D12_med <- NULL
p <- NULL
for (i in motifs) {
  print(match(i, motifs))
  motif <- c(motif, i)
  count1 <- length(which(temp$motif2 == i))
  count <- c(count, count1)
  D1.1 <- mean(temp[motif2 == i]$D1_pause)
  D1 <- c(D1, D1.1)
  D1.2 <- median(temp[motif2 == i]$D1_pause)
  D1_med <- c(D1_med, D1.2)
  D12.1 <- mean(temp[motif2 == i]$D12_pause)
  D12 <- c(D12, D12.1)
  D12.2 <- median(temp[motif2 == i]$D12_pause)
  D12_med <- c(D12_med, D12.2)
  test <- wilcox.test(temp[motif2 == i]$D1_pause, temp[motif2 == i]$D12_pause)
  p1 <- test$p.value
  p <- c(p, p1)
}
toc()
motif2_dt <- data.table(motif = motif, count = count, D1_pause = D1, D12_pause = D12,
                       D1_pause_med = D1_med, D12_pause_med = D12_med, pvalue = p)
motif2_dt[, ratio := D12_pause / D1_pause]
motif2_dt[, ratio_med := D12_pause_med / D1_pause_med]
motif2_dt[, padj := p.adjust(pvalue, method = "BH")]
saveRDS(motif2_dt, "motif2_dt.rds")
length(motif2_dt[count >= 500]$motif) # 500 used in eIF5a paper

ggplot(motif2_dt[count >= 500], aes(D1_pause, D12_pause)) + geom_point() + xlim(1.75,2.25) + ylim(1.75,2.25) + 
  geom_point(data = motif2_dt[count >= 500 & ratio > 1 & D12_pause > 2], aes(D1_pause, D12_pause), color = "blue") + 
  geom_point(data = motif2_dt[count >= 500 & ratio > 1 & D12_pause > 2 & pvalue < 0.05], aes(D1_pause, D12_pause), color = "red") + 
  geom_abline(slope = 1, intercept = 0, color = 'red') 

weblogo(motif2_dt[count >= 500 & ratio > 1 & D12_pause > 2]$motif, open = TRUE, file.out = "motif2.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
weblogo(motif2_dt[count >= 500 & ratio > 1 & D12_pause > 2]$motif, open = TRUE, file.out = "motif2_IC.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        composition = c('A'=6.325907,'C'=2.076066,'D'=5.304005,'E'=6.526124,'F'=4.787855,'G'=5.344524,'H'=2.281571,'I'=6.213785,'K'=6.339429,'L'=8.643033,'M'=2.648252,'N'=4.877476,'P'=4.892958,'Q'=4.081622,'R'=5.137781,'S'=8.080618,'T'=5.891986,'V'=6.229459,'W'=1.107328,'Y'=3.210222),
        color.scheme = 'chemistry', errorbars = FALSE, yaxis = 1)
logomaker(motif2_dt[count >= 500 & ratio > 1 & D12_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)

write.csv(motif2_dt[count >= 500 & ratio > 1 & D12_pause > 2], "/Users/KevinStein/Desktop/wormMotifs2_ratio1_D12_2.csv")
