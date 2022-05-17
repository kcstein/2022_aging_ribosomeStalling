### References ###
# http://bioconductor.org/packages/devel/bioc/vignettes/Logolas/inst/doc/Logolas.html
# https://omarwagih.github.io/ggseqlogo/
# http://kplogo.wi.mit.edu/index.html?
# https://dongyuexie.github.io/Sequence-Motif-Depletion/index.html
# https://iomics.ugent.be/icelogoserver/


library(data.table)
library(ggplot2)
devtools::install_github("collectivemedia/tictoc")
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
WT_stalls_seq <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/Chronological/doc/Fishers/WT_stalls_seq.csv", header = TRUE, stringsAsFactors = FALSE))
WT_stalls_seq60 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/Chronological/doc/Fishers/WT_stalls_60mer_seq.csv", header = TRUE, stringsAsFactors = FALSE))
sch9_stalls_seq <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/Chronological/doc/Fishers/sch9_stalls_seq.csv", header = TRUE, stringsAsFactors = FALSE))
sch9_stalls_seq60 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/Chronological/doc/Fishers/sch9_stalls_60mer_seq.csv", header = TRUE, stringsAsFactors = FALSE))
WT_stalls_seq[, tunnel := substring(sequence,1,27)]
WT_stalls_seq[, active := substring(sequence,28,30)]
WT_stalls_seq60[, nascent := substring(sequence,1,29)]
WT_stalls_seq60[, tunnel := substring(sequence,30,57)]
WT_stalls_seq60[, active := substring(sequence,58,60)]
sch9_stalls_seq[, tunnel := substring(sequence,1,27)]
sch9_stalls_seq[, active := substring(sequence,28,30)]
sch9_stalls_seq60[, nascent := substring(sequence,1,29)]
sch9_stalls_seq60[, tunnel := substring(sequence,30,57)]
sch9_stalls_seq60[, active := substring(sequence,58,60)]


### Logos ###
# for use with kplogo: A:0.05522952,C:0.01270852,D:0.05857607,E:0.06594403,F:0.04468158,G:0.05022744,H:0.02124895,I:0.06552491,K:0.07355956,L:0.09566521,M:0.02076128,N:0.06116999,P:0.04305590,Q:0.03925258,R:0.04445107,S:0.08950417,T:0.05846845,V:0.05608196,W:0.01043498,Y:0.03345382
# A=5.522952,C=1.270852,D=5.857607,E=6.594403,F=4.468158,G=5.022744,H=2.124895,I=6.552491,K=7.355956,L=9.566521,M=2.076128,N=6.116999,P=4.305590,Q=3.925258,R=4.445107,S=8.950417,T=5.846845,V=5.608196,W=1.043498,Y=3.345382
WT_stalls_seq <- readRDS("WT_stalls_seq.rds")
#bg <- c(A=0.05522952,C=0.01270852,D=0.05857607,E=0.06594403,F=0.04468158,G=0.05022744,H=0.02124895,I=0.06552491,K=0.07355956,L=0.09566521,M=0.02076128,N=0.06116999,P=0.04305590,Q=0.03925258,R=0.04445107,S=0.08950417,T=0.05846845,V=0.05608196,W=0.01043498,Y=0.03345382)
bg <- c(A=0.05487796,C=0.01280048,D=0.05896651,E=0.06613241,F=0.04491373,
        G=0.05048020,H=0.02136458,I=0.06611500,K=0.07268881,L=0.09615862,
        M=0.01905851,N=0.06164449,P=0.04329378,Q=0.03928379,R=0.04379378,
        S=0.08932823,T=0.05878230,V=0.05609427,W=0.01047119,Y=0.03375137)

weblogo(WT_stalls_seq[WT_D4_pause > 6]$active, open = TRUE, file.out = "stalling_fishers.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
weblogo(WT_stalls_seq[WT_D4_pause > 6]$active, open = TRUE, file.out = "WT_D4_prob.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        composition = c(A=5.522952,C=1.270852,D=5.857607,E=6.594403,F=4.468158,G=5.022744,H=2.124895,I=6.552491,K=7.355956,L=9.566521,M=2.076128,N=6.116999,P=4.305590,Q=3.925258,R=4.445107,S=8.950417,T=5.846845,V=5.608196,W=1.043498,Y=3.345382),
        color.scheme = 'chemistry', errorbars = FALSE, yaxis = 1)
logomaker(WT_stalls_seq$active, type = "EDLogo", bg = bg)
logomaker(WT_stalls_seq[WT_D4_pause > 10]$active, type = "EDLogo", bg = bg)
logomaker(WT_stalls_seq[WT_D4_pause > 10]$active, type = "EDLogo", bg = bg, 
          color_seed = 6)

bg <- summary(sc_fishers[position > 20 & stopdist < -20]$residue) / 
  length(sc_fishers[position > 20 & stopdist < -20]$residue)
bg <- bg[c(1:19,21)]

logomaker(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$motif3, type = "EDLogo", bg = bg, color_seed = 6)
summary(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$residue) / 
  length(WT_stalling_peaks_ageDep[WT_D4_pause > 6]$residue)

length(WT_stalling_peaks_ageDep$motif3)

bg0 <- summary(D0_Zpeaks$residue) / length(D0_Zpeaks$residue)
bg0 <- bg0[c(1:19,21)]
bg2 <- summary(D2_Zpeaks$residue) / length(D2_Zpeaks$residue)
bg2 <- bg2[c(1:19,21)]
logomaker(D2_Zpeaks$motif3, type = "EDLogo", bg = bg0, color_seed = 6)
logomaker(D4_Zpeaks$motif3, type = "EDLogo", bg = bg2, color_seed = 6)
logomaker(D4_Zpeaks$motif3, type = "EDLogo", bg = bg0, color_seed = 6)

### Genome-wide tripeptide active site analysis ###
# temp <- sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 & 
#                  WT_D2_1A_rpc >= 0.5 & WT_D2_2A_rpc >= 0.5 & 
#                  WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 & 
#                  WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 & 
#                  WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 & 
#                  WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & 
#                  WT_D0_1A > 0 & !is.na(WT_D0_1A_pause) & WT_D0_2A > 0 & !is.na(WT_D0_2A_pause) &
#                  WT_D2_1A > 0 & !is.na(WT_D2_1A_pause) & WT_D2_2A > 0 & !is.na(WT_D2_2A_pause) &
#                  WT_D4_1A > 0 & !is.na(WT_D4_1A_pause) & WT_D4_2A > 0 & !is.na(WT_D4_2A_pause) &
#                  position > 15 & stopdist < -15]
# temp[, WT_D2_pause := (WT_D2_1A_pause + WT_D2_2A_pause) / 2]
temp <- sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 & 
                 WT_D2_1A_rpc >= 0.5 & WT_D2_2A_rpc >= 0.5 & 
                 WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 & 
                 WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 & 
                 WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 & 
                 WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & 
                 WT_D0_1A > 0 & !is.na(WT_D0_1A_pause) & WT_D0_2A > 0 & !is.na(WT_D0_2A_pause) &
                 WT_D2_1A > 0 & !is.na(WT_D2_1A_pause) & WT_D2_2A > 0 & !is.na(WT_D2_2A_pause) &
                 WT_D4_1A > 0 & !is.na(WT_D4_1A_pause) & WT_D4_2A > 0 & !is.na(WT_D4_2A_pause) &
                 position > 20 & stopdist < -20]
motifs <- unique(temp[!is.na(motif3)]$motif3)

tic("runtime")
motif <- NULL
count <- NULL
D0 <- NULL
D0_med <- NULL
D2 <- NULL
D2_med <- NULL
D4 <- NULL
D4_med <- NULL
p_4_0 <- NULL
p_2_0 <- NULL
p_4_2 <- NULL
for (i in motifs) {
  print(match(i, motifs))
  motif <- c(motif, i)
  count1 <- length(which(temp$motif3 == i))
  count <- c(count, count1)
  D0.1 <- mean(temp[motif3 == i]$WT_D0_pause)
  D0 <- c(D0, D0.1)
  D0.2 <- median(temp[motif3 == i]$WT_D0_pause)
  D0_med <- c(D0_med, D0.2)
  D2.1 <- mean(temp[motif3 == i]$WT_D2_pause)
  D2 <- c(D2, D2.1)
  D2.2 <- median(temp[motif3 == i]$WT_D2_pause)
  D2_med <- c(D2_med, D2.2)
  D4.1 <- mean(temp[motif3 == i]$WT_D4_pause)
  D4 <- c(D4, D4.1)
  D4.2 <- median(temp[motif3 == i]$WT_D4_pause)
  D4_med <- c(D4_med, D4.2)
  test_4_0 <- wilcox.test(temp[motif3 == i]$WT_D0_pause, temp[motif3 == i]$WT_D4_pause)
  p4.0 <- test_4_0$p.value
  p_4_0 <- c(p_4_0, p4.0)
  test_2_0 <- wilcox.test(temp[motif3 == i]$WT_D0_pause, temp[motif3 == i]$WT_D2_pause)
  p2.0 <- test_2_0$p.value
  p_2_0 <- c(p_2_0, p2.0)
  test_4_2 <- wilcox.test(temp[motif3 == i]$WT_D2_pause, temp[motif3 == i]$WT_D4_pause)
  p4.2 <- test_4_2$p.value
  p_4_2 <- c(p_4_2, p4.2)
}
toc()
motif_dt <- data.table(motif = motif, count = count, 
                       D0_pause = D0, D2_pause = D2, D4_pause = D4,
                       D0_pause_med = D0_med, D2_pause_med = D2_med, D4_pause_med = D4_med, 
                       p_4_0 = p_4_0, p_2_0 = p_2_0, p_4_2 = p_4_2)
motif_dt[, ratio4_0 := D4_pause / D0_pause]
motif_dt[, ratio2_0 := D2_pause / D0_pause]
motif_dt[, ratio4_2 := D4_pause / D2_pause]
motif_dt[, ratio_med_4_0 := D4_pause_med / D0_pause_med]
motif_dt[, ratio_med_2_0 := D2_pause_med / D0_pause_med]
motif_dt[, ratio_med_4_2 := D4_pause_med / D2_pause_med]
motif_dt[, padj_4_0 := p.adjust(p_4_0, method = "BH")]
motif_dt[, padj_2_0 := p.adjust(p_2_0, method = "BH")]
motif_dt[, padj_4_2 := p.adjust(p_4_2, method = "BH")]
#saveRDS(motif_dt, "motif_dt.rds")
test4_0 <- motif_dt[count >= 100 & ratio4_0 > 1] # 100 used in eIF5a paper
test4_0 <- test4_0[order(D4_pause, decreasing = TRUE),]
test4_2 <- motif_dt[count >= 100 & ratio4_2 > 1] # 100 used in eIF5a paper
test4_2 <- test4_2[order(D4_pause, decreasing = TRUE),]
test2_0 <- motif_dt[count >= 100 & ratio2_0 > 1] # 100 used in eIF5a paper
test2_0 <- test2_0[order(D2_pause, decreasing = TRUE),]

tic("runtime")
motif <- NULL
count <- NULL
D0 <- NULL
D0_med <- NULL
D4 <- NULL
D4_med <- NULL
p <- NULL
for (i in motifs) {
  print(match(i, motifs))
  motif <- c(motif, i)
  count1 <- length(which(temp$motif3 == i))
  count <- c(count, count1)
  D0.1 <- mean(temp[motif3 == i]$WT_D0_pause)
  D0 <- c(D0, D0.1)
  D0.2 <- median(temp[motif3 == i]$WT_D0_pause)
  D0_med <- c(D0_med, D0.2)
  D4.1 <- mean(temp[motif3 == i]$WT_D4_pause)
  D4 <- c(D4, D4.1)
  D4.2 <- median(temp[motif3 == i]$WT_D4_pause)
  D4_med <- c(D4_med, D4.2)
  test <- wilcox.test(temp[motif3 == i]$WT_D0_pause, temp[motif3 == i]$WT_D4_pause)
  p1 <- test$p.value
  p <- c(p, p1)
}
toc()
motif_dt <- data.table(motif = motif, count = count, D0_pause = D0, D4_pause = D4,
                       D0_pause_med = D0_med, D4_pause_med = D4_med, pvalue = p)
motif_dt[, ratio := D4_pause / D0_pause]
motif_dt[, ratio_med := D4_pause_med / D0_pause_med]
motif_dt[, padj := p.adjust(pvalue, method = "BH")]
#saveRDS(motif_dt, "motif_dt.rds")

motif_dt <- readRDS("doc/motif_dt.rds")
weblogo(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 2]$motif, open = TRUE, file.out = "winner_Count_p0.05_D4above2_ratio1.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability') [D4_pause > 2]
weblogo(motif_dt[count >= 100 & ratio > 1 & D4_pause > 2]$motif, open = TRUE, file.out = "winner_Count_D4pause2_ratio1.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
# weblogo(test4_0[1:100]$motif, open = TRUE, file.out = "4_0.pdf",
#         format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
#         color.scheme = 'chemistry', units = 'probability')
# weblogo(test4_2[1:100]$motif, open = TRUE, file.out = "4_2.pdf",
#         format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
#         color.scheme = 'chemistry', units = 'probability')
# weblogo(test2_0[1:100]$motif, open = TRUE, file.out = "2_0.pdf",
#         format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
#         color.scheme = 'chemistry', units = 'probability')
weblogo(motif_dt[count >= 100 & ratio > 1 & pvalue < 0.05]$motif, open = TRUE, file.out = "Count_p_IC.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        composition = c(A=5.522952,C=1.270852,D=5.857607,E=6.594403,F=4.468158,G=5.022744,H=2.124895,I=6.552491,K=7.355956,L=9.566521,M=2.076128,N=6.116999,P=4.305590,Q=3.925258,R=4.445107,S=8.950417,T=5.846845,V=5.608196,W=1.043498,Y=3.345382),
        color.scheme = 'chemistry', errorbars = FALSE, yaxis = 1)
test1.5 <- motif_dt[count >= 100 & ratio > 1 & D4_pause > 1.5]
test1.5[, pos1 := as.factor(substring(motif,1,1))]
test1.5[, pos2 := as.factor(substring(motif,2,2))]
test1.5[, pos3 := as.factor(substring(motif,3,3))]
test2 <- motif_dt[count >= 100 & ratio > 1 & D4_pause > 2]
test2[, pos1 := as.factor(substring(motif,1,1))]
test2[, pos2 := as.factor(substring(motif,2,2))]
test2[, pos3 := as.factor(substring(motif,3,3))]
bg2 <- bg[c(1,3:6,8:18,20)]
bg1.5 <- bg[c(1,3:18,20)]
logomaker(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 2]$motif, type = "EDLogo", bg = bg, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 2]$motif, type = "EDLogo", bg = bg2, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 5)
logomaker(motif_dt[count >= 100 & ratio2_0 > 1 & D2_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 5)
logomaker(motif_dt[count >= 100 & ratio4_2 > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_2, color_seed = 5)


motif_dt[, E := as.factor(substring(motif,1,1))]
motif_dt[, P := as.factor(substring(motif,2,2))]
motif_dt[, A := as.factor(substring(motif,3,3))]
summary(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 1.5]$E)
summary(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 1.5]$P)
summary(motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 1.5]$A)
bg4_0 <- bg[c(1,3:6,8:10,12:18,20)]
summary(motif_dt[count >= 100 & ratio2_0 > 1 & D2_pause > 1.5]$E)
summary(motif_dt[count >= 100 & ratio2_0 > 1 & D2_pause > 1.5]$P)
summary(motif_dt[count >= 100 & ratio2_0 > 1 & D2_pause > 1.5]$A)
summary(motif_dt[count >= 100 & ratio4_2 > 1 & D4_pause > 1.5]$E)
summary(motif_dt[count >= 100 & ratio4_2 > 1 & D4_pause > 1.5]$P)
summary(motif_dt[count >= 100 & ratio4_2 > 1 & D4_pause > 1.5]$A)
bg4_2 <- bg[c(1,3:6,8:18,20)]



motif4_0 <- motif_dt[count >= 100 & ratio4_0 > 1 & D4_pause > 1.5]
motif2_0 <- motif_dt[count >= 100 & ratio2_0 > 1 & D2_pause > 1.5]
motif4_2 <- summary(motif_dt[count >= 100 & ratio4_2 > 1 & D4_pause > 1.5]$A)

D2freq <- (summary(motif2_0$A) / length(motif2_0$A)) / bg
  (summary(motif_dt[count >= 100]$A) / length(motif_dt[count >= 100]$A))
D4freq <- (summary(motif4_2$A) / length(motif4_2$A)) / bg
  (summary(motif_dt[count >= 100]$A) / length(motif_dt[count >= 100]$A))
freq2 <- data.table(age = "D2", residue = names(D2freq), freq = D2freq)
freq4 <- data.table(age = "D4", residue = names(D4freq), freq = D4freq)
freq <- rbind(freq2, freq4)
freq$age <- factor(freq$age, levels = c("D2", "D4"))

ggplot(freq[residue != "X"], aes(residue,log2(freq), fill = age)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(limits = c("D2", "D4"), 
                    labels = c("Day 2", "Day 4"),
                    values = c("#999999", "#E31A1C"), name = "")


# to use with kplogo online server
write.csv(motif_dt[count >= 100 & ratio > 1 & D4_pause > 1.5], "/Users/KevinStein/Desktop/yeast_motifs1.5.csv")


ggplot(motif_dt[count >= 100], aes(D0_pause, D4_pause)) + 
  geom_abline(intercept = 0, slope = 1, color = "blue") + geom_point() +
  geom_text(data = motif_dt[count >= 100 & ratio4_0 > 1.4 & D4_pause > 1.5], aes(D0_pause, D4_pause, label = motif), color = "red")

ggplot(collisions[adjusted == 0], aes(WT_D0_pause, WT_D4_pause)) + 
  geom_abline(intercept = 0, slope = 1, color = "blue") + geom_point() +
  geom_text(data = collisions[adjusted == 0], aes(WT_D0_pause, WT_D4_pause, label = motif3), color = "red")

### Genome-wide tripeptide active site analysis ###
temp <- sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 & 
                 WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 & 
                 WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 & 
                 WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & 
                 WT_D0_1A > 0 & !is.na(WT_D0_1A_pause) & WT_D0_2A > 0 & !is.na(WT_D0_2A_pause) &
                 WT_D4_1A > 0 & !is.na(WT_D4_1A_pause) & WT_D4_2A > 0 & !is.na(WT_D4_2A_pause) &
                 position > 15 & stopdist < -15]
motifs <- unique(temp[!is.na(motif2)]$motif2)

tic("runtime")
motif <- NULL
count <- NULL
D0 <- NULL
D0_med <- NULL
D4 <- NULL
D4_med <- NULL
p <- NULL
for (i in motifs) {
  print(match(i, motifs))
  motif <- c(motif, i)
  count1 <- length(which(temp$motif2 == i))
  count <- c(count, count1)
  D0.1 <- mean(temp[motif2 == i]$WT_D0_pause)
  D0 <- c(D0, D0.1)
  D0.2 <- median(temp[motif2 == i]$WT_D0_pause)
  D0_med <- c(D0_med, D0.2)
  D4.1 <- mean(temp[motif2 == i]$WT_D4_pause)
  D4 <- c(D4, D4.1)
  D4.2 <- median(temp[motif2 == i]$WT_D4_pause)
  D4_med <- c(D4_med, D4.2)
  test <- wilcox.test(temp[motif2 == i]$WT_D0_pause, temp[motif2 == i]$WT_D4_pause)
  p1 <- test$p.value
  p <- c(p, p1)
}
toc()
motif2_dt <- data.table(motif = motif, count = count, D0_pause = D0, D4_pause = D4,
                       D0_pause_med = D0_med, D4_pause_med = D4_med, pvalue = p)
motif2_dt[, ratio := D4_pause / D0_pause]
motif2_dt[, ratio_med := D4_pause_med / D0_pause_med]
motif2_dt[, padj := p.adjust(pvalue, method = "BH")]
saveRDS(motif2_dt, "motif2_dt.rds")

weblogo(motif2_dt[count >= 500 & ratio > 1 & D4_pause > 1.5]$motif, open = TRUE, file.out = "dipeptide_prob.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
weblogo(motif2_dt[count >= 500 & ratio > 1 & D4_pause > 1.5]$motif, open = TRUE, file.out = "dipeptide_IC.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        composition = c(A=5.522952,C=1.270852,D=5.857607,E=6.594403,F=4.468158,G=5.022744,H=2.124895,I=6.552491,K=7.355956,L=9.566521,M=2.076128,N=6.116999,P=4.305590,Q=3.925258,R=4.445107,S=8.950417,T=5.846845,V=5.608196,W=1.043498,Y=3.345382),
        color.scheme = 'chemistry', errorbars = FALSE, yaxis = 1.6)

# to use with kplogo online server
write.csv(motif2_dt[count >= 500 & ratio > 1 & D4_pause > 1.5], "/Users/KevinStein/Desktop/yeast_dipeptide.csv")



### Compare Day 4 WT and sch9
temp <- sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 & 
                 WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 & 
                 WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 & 
                 WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & 
                 WT_D0_1A > 0 & !is.na(WT_D0_1A_pause) & WT_D0_2A > 0 & !is.na(WT_D0_2A_pause) &
                 WT_D4_1A > 0 & !is.na(WT_D4_1A_pause) & WT_D4_2A > 0 & !is.na(WT_D4_2A_pause) &
                 sch9_D0_1A_rpc >= 0.5 & sch9_D0_2A_rpc >= 0.5 & 
                 sch9_D4_1A_rpc >= 0.5 & sch9_D4_2A_rpc >= 0.5 & 
                 sch9_D0_1A_sum >= 64 & sch9_D0_2A_sum >= 64 & 
                 sch9_D4_1A_sum >= 64 & sch9_D4_2A_sum >= 64 & 
                 sch9_D0_1A > 0 & !is.na(sch9_D0_1A_pause) & sch9_D0_2A > 0 & !is.na(sch9_D0_2A_pause) &
                 sch9_D4_1A > 0 & !is.na(sch9_D4_1A_pause) & sch9_D4_2A > 0 & !is.na(sch9_D4_2A_pause) &
                 position > 20 & stopdist < -20]
motifs <- unique(temp[!is.na(motif3)]$motif3)

tic("runtime")
motif <- NULL
count <- NULL
sch9_D4 <- NULL
D4 <- NULL
for (i in motifs) {
  print(match(i, motifs))
  motif <- c(motif, i)
  count1 <- length(which(temp$motif3 == i))
  count <- c(count, count1)
  sch9_D4.1 <- mean(temp[motif3 == i]$sch9_D4_pause)
  sch9_D4 <- c(sch9_D4, sch9_D4.1)
  D4.1 <- mean(temp[motif3 == i]$WT_D4_pause)
  D4 <- c(D4, D4.1)
}
toc()
motif_dt <- data.table(motif = motif, count = count, 
                       sch9_D4_pause = sch9_D4, D4_pause = D4)
motif_dt[, ratio := D4_pause / sch9_D4_pause]

motif_dt[, E := as.factor(substring(motif,1,1))]
motif_dt[, P := as.factor(substring(motif,2,2))]
motif_dt[, A := as.factor(substring(motif,3,3))]
summary(motif_dt[count >= 100 & ratio > 1 & D4_pause > 2]$E)
summary(motif_dt[count >= 100 & ratio > 1 & D4_pause > 2]$P)
summary(motif_dt[count >= 100 & ratio > 1 & D4_pause > 2]$A)
bg4_0 <- bg[c(1,3:18,20)]
bg2 <- bg[c(1,3:6,8:10,12:18,20)]

logomaker(motif_dt[count >= 100 & ratio > 1 & D4_pause > 1.5]$motif, type = "EDLogo", bg = bg4_0, color_seed = 6)
logomaker(motif_dt[count >= 100 & ratio > 1 & D4_pause > 2]$motif, type = "EDLogo", bg = bg2, color_seed = 6)

weblogo(motif_dt[count >= 100 & ratio > 1.3]$motif, open = TRUE, file.out = "stalling_fishers.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        color.scheme = 'chemistry', units = 'probability')
weblogo(motif_dt[count >= 100 & ratio > 1]$motif, open = TRUE, file.out = "WT_D4_prob.pdf",
        format = 'pdf', sequence.type = 'protein', alphabet = 'ACDEFGHIKLMNPQRSTVWY',
        composition = c(A=5.522952,C=1.270852,D=5.857607,E=6.594403,F=4.468158,G=5.022744,H=2.124895,I=6.552491,K=7.355956,L=9.566521,M=2.076128,N=6.116999,P=4.305590,Q=3.925258,R=4.445107,S=8.950417,T=5.846845,V=5.608196,W=1.043498,Y=3.345382),
        color.scheme = 'chemistry', errorbars = FALSE, yaxis = 1)

