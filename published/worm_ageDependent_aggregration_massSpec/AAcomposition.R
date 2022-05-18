library(data.table)
library(plyr)
library(dplyr)
library(tidyr)


### Composition ###
ce_comp <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/proteome_AAcomposition_fraction.csv', header = TRUE, stringsAsFactors = TRUE))

dir <- "/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/composition/"
files <- list.files(dir, pattern = "*.csv")
samples <- sub("\\.csv", '', files)

for(i in 1:length(files)) {
  print(files[i])
  j <- as.data.table(lapply(paste0(dir, files[i]), read.csv, header = FALSE, stringsAsFactors = TRUE))
  setnames(j, c('orf', 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
  setkeyv(j, c("orf"))
  assign(samples[i], j)
}



wilcox.test(kenyon_aggregates_AAcomposition$P, ce_comp$P, alternative = 'l')
summary(ce_comp$R)
summary(hartl_proteome_AAcomposition$R)
summary(hartl_aggregates_AAcomposition$R)
summary(hartl_aggregatesANDstall_AAcomposition$R)
summary(kenyon_aggregates_AAcomposition$R)
summary(kenyon_aggregatesANDstall_AAcomposition$R)
summary(stalling_ageDep_Proteins_AAcomposition$R)
summary(stalling_peaks_hartlProteome_Proteins_AAcomposition$R)

wilcox.test(hartl_aggregates_AAcomposition$R, ce_comp$R, alternative = 't')
wilcox.test(hartl_aggregatesANDstall_AAcomposition$R, ce_comp$R, alternative = 't')
wilcox.test(kenyon_aggregates_AAcomposition$R, ce_comp$R, alternative = 't')
wilcox.test(kenyon_aggregatesANDstall_AAcomposition$R, ce_comp$R, alternative = 't')
wilcox.test(stalling_ageDep_Proteins_AAcomposition$R, ce_comp$R, alternative = 't')
wilcox.test(stalling_peaks_hartlProteome_Proteins_AAcomposition$R, ce_comp$R, alternative = 't')

summary(N2_dtA[position > 0 & stopdist < 0]$residue) / length(N2_dtA[position > 0 & stopdist < 0]$residue) 

aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
stalled <- NULL
stalled_P10 <- NULL
hartl <- NULL
hartl_stalled <- NULL
kenyon <- NULL
kenyon_stalled <- NULL
proteome <- NULL
hartl_MSbg <- NULL
p_stalled <- NULL
p_stalled_P10 <- NULL
p_hartl <- NULL
p_hartl_stalled <- NULL
p_hartl_MSbg <- NULL
p_kenyon <- NULL
p_kenyon_stalled <- NULL
p_kenyon_MSbg <- NULL
for (i in aa) {
  stalled1 <- mean(stallingProteins[[i]])
  stalled <- c(stalled, stalled1)
  stalled_P101 <- mean(stallingProteins_P10[[i]])
  stalled_P10 <- c(stalled_P10, stalled_P101)
  hartl1 <- mean(hartl_aggregates_AAcomposition[[i]])
  hartl <- c(hartl, hartl1)
  hartl_stalled1 <- mean(hartl_aggregatesANDstalling_AAcomposition[[i]])
  hartl_stalled <- c(hartl_stalled, hartl_stalled1)
  kenyon1 <- mean(kenyon_aggregates_AAcomposition[[i]])
  kenyon <- c(kenyon, kenyon1)
  kenyon_stalled1 <- mean(kenyon_aggregatesANDstalling_AAcomposition[[i]])
  kenyon_stalled <- c(kenyon_stalled, kenyon_stalled1)
  proteome1 <- mean(ce_comp[[i]])
  proteome <- c(proteome, proteome1)
  hartl_MSbg1 <- mean(background_hartl_composition[[i]])
  hartl_MSbg <- c(hartl_MSbg, hartl_MSbg1)
  p_stalled1 <- wilcox.test(stallingProteins[[i]], 
                           ce_comp[[i]], alternative = 't')$p.value
  p_stalled <- c(p_stalled, p_stalled1)
  p_stalled_P101 <- wilcox.test(stallingProteins_P10[[i]], 
                            ce_comp[[i]], alternative = 't')$p.value
  p_stalled_P10 <- c(p_stalled_P10, p_stalled_P101)
  p_hartl1 <- wilcox.test(hartl_aggregates_AAcomposition[[i]], 
                            ce_comp[[i]], alternative = 't')$p.value
  p_hartl <- c(p_hartl, p_hartl1)
  p_hartl_stalled1 <- wilcox.test(hartl_aggregatesANDstalling_AAcomposition[[i]], 
                            ce_comp[[i]], alternative = 't')$p.value
  p_hartl_stalled <- c(p_hartl_stalled, p_hartl_stalled1)
  p_hartl_MSbg1 <- wilcox.test(hartl_aggregates_AAcomposition[[i]], 
                          background_hartl_composition[[i]], alternative = 't')$p.value
  p_hartl_MSbg <- c(p_hartl_MSbg, p_hartl_MSbg1)
  p_kenyon1 <- wilcox.test(kenyon_aggregates_AAcomposition[[i]], 
                           ce_comp[[i]], alternative = 't')$p.value
  p_kenyon <- c(p_kenyon, p_kenyon1)
  p_kenyon_stalled1 <- wilcox.test(kenyon_aggregatesANDstalling_AAcomposition[[i]], 
                                   ce_comp[[i]], alternative = 't')$p.value
  p_kenyon_stalled <- c(p_kenyon_stalled, p_kenyon_stalled1)
  p_kenyon_MSbg1 <- wilcox.test(kenyon_aggregates_AAcomposition[[i]], 
                               background_hartl_composition[[i]], alternative = 't')$p.value
  p_kenyon_MSbg <- c(p_kenyon_MSbg, p_kenyon_MSbg1)
}
composition <- data.table(aa, proteome, stalled, p_stalled, stalled_P10, p_stalled_P10, 
                          hartl_MSbg, hartl, p_hartl, p_hartl_MSbg, hartl_stalled, p_hartl_stalled,
                          kenyon, p_kenyon, p_kenyon_MSbg, kenyon_stalled, p_kenyon_stalled)
composition[, stalled_ratio := stalled / proteome]
composition[, stalled_P10_ratio := stalled_P10 / proteome]
composition[, hartl_ratio := hartl / proteome]
composition[, hartl_MSbg_ratio := hartl / hartl_MSbg]
composition[, hartl_stalled_ratio := hartl_stalled / proteome]
composition[, kenyon_ratio := kenyon / proteome]
composition[, kenyon_MSbg_ratio := kenyon / hartl_MSbg]
composition[, kenyon_stalled_ratio := kenyon_stalled / proteome]
View(composition[aa == 'R'])

ggplot(composition, aes(x=aa, y=log2(stalled_ratio))) + geom_col()
ggplot(composition, aes(x=aa, y=log2(stalled_P10_ratio))) + geom_col()
ggplot(composition, aes(x=aa, y=log2(hartl_ratio))) + geom_col()
ggplot(composition, aes(x=aa, y=log2(hartl_MSbg_ratio))) + geom_col()
ggplot(composition, aes(x=aa, y=log2(hartl_stalled_ratio))) + geom_col()
ggplot(composition, aes(x=aa, y=log2(kenyon_ratio))) + geom_col()
ggplot(composition, aes(x=aa, y=log2(kenyon_MSbg_ratio))) + geom_col()
ggplot(composition, aes(x=aa, y=log2(kenyon_stalled_ratio))) + geom_col()


stalled <- NULL
stalled_P10 <- NULL
hartl <- NULL
hartl_stalled <- NULL
kenyon <- NULL
kenyon_stalled <- NULL
proteome <- NULL
hartl_MSbg <- NULL
p_stalled <- NULL
p_stalled_P10 <- NULL
p_hartl <- NULL
p_hartl_stalled <- NULL
p_hartl_MSbg <- NULL
p_kenyon <- NULL
p_kenyon_stalled <- NULL
p_kenyon_MSbg <- NULL
for (i in aa) {
  stalled1 <- median(stallingProteins[[i]])
  stalled <- c(stalled, stalled1)
  stalled_P101 <- median(stallingProteins_P10[[i]])
  stalled_P10 <- c(stalled_P10, stalled_P101)
  hartl1 <- median(hartl_aggregates_AAcomposition[[i]])
  hartl <- c(hartl, hartl1)
  hartl_stalled1 <- median(hartl_aggregatesANDstalling_AAcomposition[[i]])
  hartl_stalled <- c(hartl_stalled, hartl_stalled1)
  kenyon1 <- median(kenyon_aggregates_AAcomposition[[i]])
  kenyon <- c(kenyon, kenyon1)
  kenyon_stalled1 <- median(kenyon_aggregatesANDstalling_AAcomposition[[i]])
  kenyon_stalled <- c(kenyon_stalled, kenyon_stalled1)
  proteome1 <- median(ce_comp[[i]])
  proteome <- c(proteome, proteome1)
  hartl_MSbg1 <- median(background_hartl_composition[[i]])
  hartl_MSbg <- c(hartl_MSbg, hartl_MSbg1)
  p_stalled1 <- wilcox.test(stallingProteins[[i]], 
                            ce_comp[[i]], alternative = 't')$p.value
  p_stalled <- c(p_stalled, p_stalled1)
  p_stalled_P101 <- wilcox.test(stallingProteins_P10[[i]], 
                                ce_comp[[i]], alternative = 't')$p.value
  p_stalled_P10 <- c(p_stalled_P10, p_stalled_P101)
  p_hartl1 <- wilcox.test(hartl_aggregates_AAcomposition[[i]], 
                          ce_comp[[i]], alternative = 't')$p.value
  p_hartl <- c(p_hartl, p_hartl1)
  p_hartl_stalled1 <- wilcox.test(hartl_aggregatesANDstalling_AAcomposition[[i]], 
                                  ce_comp[[i]], alternative = 't')$p.value
  p_hartl_stalled <- c(p_hartl_stalled, p_hartl_stalled1)
  p_hartl_MSbg1 <- wilcox.test(hartl_aggregates_AAcomposition[[i]], 
                               background_hartl_composition[[i]], alternative = 't')$p.value
  p_hartl_MSbg <- c(p_hartl_MSbg, p_hartl_MSbg1)
  p_kenyon1 <- wilcox.test(kenyon_aggregates_AAcomposition[[i]], 
                           ce_comp[[i]], alternative = 't')$p.value
  p_kenyon <- c(p_kenyon, p_kenyon1)
  p_kenyon_stalled1 <- wilcox.test(kenyon_aggregatesANDstalling_AAcomposition[[i]], 
                                   ce_comp[[i]], alternative = 't')$p.value
  p_kenyon_stalled <- c(p_kenyon_stalled, p_kenyon_stalled1)
  p_kenyon_MSbg1 <- wilcox.test(kenyon_aggregates_AAcomposition[[i]], 
                                background_hartl_composition[[i]], alternative = 't')$p.value
  p_kenyon_MSbg <- c(p_kenyon_MSbg, p_kenyon_MSbg1)
}
composition_median <- data.table(aa, proteome, stalled, p_stalled, stalled_P10, p_stalled_P10, 
                                 hartl_MSbg, hartl, p_hartl, p_hartl_MSbg, hartl_stalled, p_hartl_stalled,
                                 kenyon, p_kenyon, p_kenyon_MSbg, kenyon_stalled, p_kenyon_stalled)
composition_median[, stalled_ratio := stalled / proteome]
composition_median[, stalled_P10_ratio := stalled_P10 / proteome]
composition_median[, hartl_ratio := hartl / proteome]
composition_median[, hartl_MSbg_ratio := hartl / hartl_MSbg]
composition_median[, hartl_stalled_ratio := hartl_stalled / proteome]
composition_median[, kenyon_ratio := kenyon / proteome]
composition_median[, kenyon_MSbg_ratio := kenyon / hartl_MSbg]
composition_median[, kenyon_stalled_ratio := kenyon_stalled / proteome]
View(composition_median[aa == 'R'])
View(composition_median[stalled_ratio > 1 & stalled_P10_ratio > 1 &
                          hartl_ratio > 1 & hartl_MSbg_ratio > 1 &
                          kenyon_ratio > 1 & kenyon_MSbg_ratio > 1])

ggplot(composition_median, aes(x=aa, y=log2(stalled_ratio))) + geom_col()
ggplot(composition_median, aes(x=aa, y=log2(stalled_P10_ratio))) + geom_col()
ggplot(composition_median, aes(x=aa, y=log2(hartl_ratio))) + geom_col()
ggplot(composition_median, aes(x=aa, y=log2(hartl_MSbg_ratio))) + geom_col()
ggplot(composition_median, aes(x=aa, y=log2(hartl_stalled_ratio))) + geom_col()
ggplot(composition_median, aes(x=aa, y=log2(kenyon_ratio))) + geom_col()
ggplot(composition_median, aes(x=aa, y=log2(kenyon_MSbg_ratio))) + geom_col()
ggplot(composition_median, aes(x=aa, y=log2(kenyon_stalled_ratio))) + geom_col()