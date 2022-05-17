library(data.table)
library(dplyr)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')

# Import codon density files
dir <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Young_Green_2015/"
files <- list.files(dir, pattern = "*.codon")
samples <- sub("\\.codon", '', files)

for(i in 1:length(files)) {
  print(files[i])
  assign(samples[i], as.data.frame(lapply(paste0(dir, files[i]), read.delim, header = FALSE, stringsAsFactors = TRUE)))
}
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

# Create data table with basic ORF properties
his_dt <- data.table(orf = dfs[[1]][, 1],
                    position = dfs[[1]][, 2],
                    codon = dfs[[1]][, 3]) 
setkeyv(his_dt, c("orf"))
orfs <- unique(his_dt[, 1])
his_dt[, length := (length(position) - 15), by = orf] # subtract flanking 7 codons and stop codon
his_dt[, stopdist := position - (length + 1)] # stop codon is at stopdist == 0
his_dt[, ID := as.character(paste(orf, position, sep = "_"))]

# Add residue
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(his_dt$codon, codon_table$codon))
his_dt <- cbind(his_dt, residue = codon_table[i]$residue)

# Add raw counts
for(i in 1:length(samples)) {
  print(samples[i])
  his_dt <- cbind(his_dt, j = dfs[[i]][,4])
  setnames(his_dt, "j", samples[i])
}


# Calculate rpc for each orf after excluding first and last 5 sense codons (Rachel Green excludes 100nt in the eIF5A paper)
# Sum reads in order to exclude genes with less than 64 reads in each replicate in downstream analysis
temp <- his_dt[position > 30 & stopdist < -30]
for(i in samples) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  temp[, med1 := lapply(.SD[, ..i], median), by = orf]
  gene <- cbind(match(his_dt$orf, temp$orf))
  his_dt <- cbind(his_dt, sum1 = temp[gene]$sum1)
  his_dt <- cbind(his_dt, rpc1 = temp[gene]$rpc1)
  his_dt <- cbind(his_dt, med1 = temp[gene]$med1)
  setnames(his_dt, "sum1", paste0(i,"_sum"))
  setnames(his_dt, "rpc1", paste0(i,"_rpc"))
  setnames(his_dt, "med1", paste0(i,"_med"))
}
saveRDS(his_dt, "his_dt.rds")


# Calculate rpm, pause score, and z-score
#his_dt <- readRDS("his_dt.rds")
for(i in samples) {
  print(i)
  his_dt[, paste0(i,"_rpm") := (his_dt[[i]] / sum(his_dt[[i]])) * 10^6] # calculate rpm
  his_dt[, paste0(i,"_pause") := his_dt[[i]] / his_dt[[paste0(i,"_rpc")]]] # calculate pause
  his_dt[, Z := abs(his_dt[[i]] - his_dt[[paste0(i,"_med")]])]
  his_dt[, Z := median(Z), by = orf]
  his_dt[, Z := (0.6745*(his_dt[[i]] - his_dt[[paste0(i,"_med")]])) / Z] # calculate z-score
  setnames(his_dt, "Z", paste0(i,"_Zmad"))
  his_dt[, paste0(i,"_Z") := (his_dt[[i]] - his_dt[[paste0(i,"_rpc")]]) / ((his_dt[[paste0(i,"_rpc")]])^(1/2))]
}
saveRDS(his_dt, "his_dt.rds")


### Test
his_dt <- readRDS("his_dt.rds")
ggplot(his_dt[WT1_A_rpc >= 1 & AT1_A_rpc >= 1 & residue != "X"], aes(residue, WT1_A_Z)) + geom_violin(aes(residue, AT1_A_Z), color = "red", fill = NA) + geom_violin(fill = NA) + 
  scale_y_log10(limits = c(0.01,100)) # doesn't make sense to plot Z scores since half under 0
ggplot(his_dt[WT2_A_rpc >= 1 & AT2_A_rpc >= 1 & residue != "X"], aes(residue, WT2_A_Z)) + geom_violin(aes(residue, AT2_A_Z), color = "red", fill = NA) + geom_violin(fill = NA) + 
  scale_y_log10(limits = c(0.01,100))
ggplot(his_dt[WT1_A_rpc >= 1 & AT1_A_rpc >= 1 & residue != "X"], aes(residue, WT1_A_pause)) + geom_violin(aes(residue, AT1_A_pause), color = "red", fill = NA) + geom_violin(fill = NA) +
  scale_y_log10(limits = c(0.01,100))
ggplot(his_dt[WT1_A_rpc >= 1 & AT1_A_rpc >= 1 & residue != "X"], aes(residue, WT1_P_pause)) + geom_violin(aes(residue, AT1_P_pause), color = "red", fill = NA) + geom_violin(fill = NA) +
  scale_y_log10(limits = c(0.01,100))
ggplot(his_dt[WT1_A_rpc >= 1 & AT1_A_rpc >= 1 & residue != "X"], aes(residue, WT1_E_pause)) + geom_violin(aes(residue, AT1_E_pause), color = "red", fill = NA) + geom_violin(fill = NA) +
  scale_y_log10(limits = c(0.01,100))
ggplot(his_dt[WT2_A_rpc >= 1 & AT2_A_rpc >= 1 & residue != "X"], aes(residue, WT2_A_pause)) + geom_violin(aes(residue, AT2_A_pause), color = "red", fill = NA) + geom_violin(fill = NA) +
  scale_y_log10(limits = c(0.01,100))

ggplot(his_dt[WT2_A_rpc >= 0.5 & AT2_A_rpc >= 0.5], aes(WT2_A_pause)) + stat_ecdf(geom = "step", color = "black") +
  scale_x_log10(limits = c(0.01,100)) +
  stat_ecdf(aes(AT2_A_pause), geom = "step", color = "red")

AT_Ponly <- AT_P10[!AT_P10$ID %in% ATzmad10A$ID]
AT_Zonly <- ATzmad10A[!ATzmad10A$ID %in% AT_P10$ID]
AT_combined <- AT_P10[AT_P10$ID %in% ATzmad10A$ID]

ggplot(AT_P10, aes(AT2_A_pause)) + stat_ecdf(geom = "step", color = "red") +
  scale_x_log10(limits = c(0.01,100)) +
  stat_ecdf(aes(AT1_A_pause), geom = "step", color = "red") +
  stat_ecdf(data = ATzmad10A, aes(AT1_A_pause), geom = "step", color = "black") +
  stat_ecdf(data = ATzmad10A, aes(AT2_A_pause), geom = "step", color = "black") +
  stat_ecdf(data = AT_Ponly, aes(AT1_A_pause), geom = "step", color = "orange") +
  stat_ecdf(data = AT_Ponly, aes(AT2_A_pause), geom = "step", color = "orange") +
  stat_ecdf(data = AT_Zonly, aes(AT1_A_pause), geom = "step", color = "blue") +
  stat_ecdf(data = AT_Zonly, aes(AT2_A_pause), geom = "step", color = "blue") +
  stat_ecdf(data = AT_combined, aes(AT1_A_pause), geom = "step", color = "green") +
  stat_ecdf(data = AT_combined, aes(AT2_A_pause), geom = "step", color = "green")
# The plot above demonstrates that sites identified using Zmad include sites that have
# a low pause score. This is likely due to these sites having statistically
# significant greater coverage. By corallary, sites identified by pause score and not Zmad
# likely have greater variation in ribosome density.


WTz10A <- his_dt[WT1_A_rpc >= 0.5 & WT2_A_rpc >= 0.5 & residue != "X" & position > 30 & WT1_A_Z >= 10 & WT2_A_Z >= 10]
WTzmad10A <- his_dt[WT1_A_rpc >= 0.5 & WT2_A_rpc >= 0.5 & residue != "X" & position > 30 & WT1_A_Zmad >= 10 & WT2_A_Zmad >= 10 & WT1_A_Zmad < Inf & WT2_A_Zmad < Inf]
WTz3A <- his_dt[WT1_A_rpc >= 0.5 & WT2_A_rpc >= 0.5 & residue != "X" & position > 30 & WT1_A_Z > 3.5 & WT2_A_Z > 3.5]
WTzmad3A <- his_dt[WT1_A_rpc >= 0.5 & WT2_A_rpc >= 0.5 & residue != "X" & position > 30 & WT1_A_Zmad > 3.5 & WT2_A_Zmad > 3.5 & WT1_A_Zmad < Inf & WT2_A_Zmad < Inf]
ATz10A <- his_dt[AT1_A_rpc >= 0.5 & AT2_A_rpc >= 0.5 & residue != "X" & position > 30 & AT1_A_Z >= 10 & AT2_A_Z >= 10]
ATzmad10A <- his_dt[AT1_A_rpc >= 0.5 & AT2_A_rpc >= 0.5 & residue != "X" & position > 30 & AT1_A_Zmad >= 10 & AT2_A_Zmad >= 10 & AT1_A_Zmad < Inf & AT2_A_Zmad < Inf]
ATz3A <- his_dt[AT1_A_rpc >= 0.5 & AT2_A_rpc >= 0.5 & residue != "X" & position > 30 & AT1_A_Z > 3.5 & AT2_A_Z > 3.5]
ATzmad3A <- his_dt[AT1_A_rpc >= 0.5 & AT2_A_rpc >= 0.5 & residue != "X" & position > 30 & AT1_A_Zmad > 3.5 & AT2_A_Zmad > 3.5 & AT1_A_Zmad < Inf & AT2_A_Zmad < Inf]
WTz3E <- his_dt[WT1_E_rpc >= 0.5 & WT2_E_rpc >= 0.5 & residue != "X" & position > 30 & WT1_E_Z > 3.5 & WT2_E_Z > 3.5]
ATz3E <- his_dt[AT1_E_rpc >= 0.5 & AT2_E_rpc >= 0.5 & residue != "X" & position > 30 & AT1_E_Z > 3.5 & AT2_E_Z > 3.5]
WTz3P <- his_dt[WT1_P_rpc >= 0.5 & WT2_P_rpc >= 0.5 & residue != "X" & position > 30 & WT1_P_Z > 3.5 & WT2_P_Z > 3.5]
ATz3P <- his_dt[AT1_P_rpc >= 0.5 & AT2_P_rpc >= 0.5 & residue != "X" & position > 30 & AT1_P_Z > 3.5 & AT2_P_Z > 3.5]
WT_P10 <- his_dt[WT1_A_rpc >= 0.5 & WT2_A_rpc >= 0.5 & residue != "X" & position > 30 & WT1_A_pause >= 10 & WT2_A_pause >= 10]
AT_P10 <- his_dt[AT1_A_rpc >= 0.5 & AT2_A_rpc >= 0.5 & residue != "X" & position > 30 & AT1_A_pause >= 10 & AT2_A_pause >= 10]

# Comparing A-site vs P and E
residue_freq <- data.table(WT10A = table(WTz10A$residue),
                           AT10A = table(ATz10A$residue),
                           WT3A = table(WTz3A$residue),
                           AT3A = table(ATz3A$residue),
                           WT3P = table(WTz3P$residue),
                           AT3P = table(ATz3P$residue),
                           WT3E = table(WTz3E$residue),
                           AT3E = table(ATz3E$residue))
residue_freq <- data.table(residue = residue_freq[, 1],
                           WT10A = (residue_freq[, 2] + 1),
                           AT10A = (residue_freq[, 4] + 1),
                           WT3A = (residue_freq[, 6] + 1),
                           AT3A = (residue_freq[, 8] + 1),
                           WT3P = (residue_freq[, 10] + 1),
                           AT3P = (residue_freq[, 12] + 1),
                           WT3E = (residue_freq[, 14] + 1),
                           AT3E = (residue_freq[, 16] + 1))
setnames(residue_freq, c("residue", "WT10A", "AT10A", "WT3A", "AT3A", "WT3P", "AT3P", "WT3E", "AT3E"))

residue_freq[, WT10A_freq := (WT10A / sum(WT10A))*100]
residue_freq[, AT10A_freq := (AT10A / sum(AT10A))*100]
residue_freq[, WT3A_freq := (WT3A / sum(WT3A))*100]
residue_freq[, AT3A_freq := (AT3A / sum(AT3A))*100]
residue_freq[, WT3P_freq := (WT3P / sum(WT3P))*100]
residue_freq[, AT3P_freq := (AT3P / sum(AT3P))*100]
residue_freq[, WT3E_freq := (WT3E / sum(WT3E))*100]
residue_freq[, AT3E_freq := (AT3E / sum(AT3E))*100]

residue_freq[, diff10A := AT10A_freq - WT10A_freq]
residue_freq[, diff3A := AT3A_freq - WT3A_freq]
residue_freq[, div10A := AT10A_freq / WT10A_freq]
residue_freq[, div3A := AT3A_freq / WT3A_freq]
residue_freq[, div3P := AT3P_freq / WT3P_freq]
residue_freq[, div3E := AT3E_freq / WT3E_freq]

plot <- ggplot(residue_freq, aes(residue, diff10A)) + geom_col()
ggplot(residue_freq, aes(residue, diff3A)) + geom_col()
ggplot(residue_freq, aes(residue, log2(div10A))) + geom_col()
ggplot(residue_freq, aes(residue, log2(div3A))) + geom_col()
ggplot(residue_freq, aes(residue, log2(div3P))) + geom_col()
ggplot(residue_freq, aes(residue, log2(div3E))) + geom_col()


# Comparing Z scores
residue_freq <- data.table(WT10Z = table(WTz10A$residue),
                           WT10Zmad = table(WTzmad10A$residue),
                           WT3Z = table(WTz3A$residue),
                           WT3Zmad = table(WTzmad3A$residue),
                           AT10Z = table(ATz10A$residue),
                           AT10Zmad = table(ATzmad10A$residue),
                           AT3Z = table(ATz3A$residue),
                           AT3Zmad = table(ATzmad3A$residue))
residue_freq <- data.table(residue = residue_freq[, 1],
                           WT10Z = (residue_freq[, 2] + 1),
                           WT10Zmad = (residue_freq[, 4] + 1),
                           WT3Z = (residue_freq[, 6] + 1),
                           WT3Zmad = (residue_freq[, 8] + 1),
                           AT10Z = (residue_freq[, 10] + 1),
                           AT10Zmad = (residue_freq[, 12] + 1),
                           AT3Z = (residue_freq[, 14] + 1),
                           AT3Zmad = (residue_freq[, 16] + 1))
setnames(residue_freq, c("residue", "WT10Z", "WT10Zmad", "WT3Z", "WT3Zmad",
                         "AT10Z", "AT10Zmad", "AT3Z", "AT3Zmad"))

residue_freq[, WT10Z_freq := (WT10Z / sum(WT10Z))*100]
residue_freq[, WT10Zmad_freq := (WT10Zmad / sum(WT10Zmad))*100]
residue_freq[, WT3Z_freq := (WT3Z / sum(WT3Z))*100]
residue_freq[, WT3Zmad_freq := (WT3Zmad / sum(WT3Zmad))*100]
residue_freq[, AT10Z_freq := (AT10Z / sum(AT10Z))*100]
residue_freq[, AT10Zmad_freq := (AT10Zmad / sum(AT10Zmad))*100]
residue_freq[, AT3Z_freq := (AT3Z / sum(AT3Z))*100]
residue_freq[, AT3Zmad_freq := (AT3Zmad / sum(AT3Zmad))*100]

residue_freq[, div10Z := AT10Z_freq / WT10Z_freq]
residue_freq[, div10Zmad := AT10Zmad_freq / WT10Zmad_freq]
residue_freq[, div3Z := AT3Z_freq / WT3Z_freq]
residue_freq[, div3Zmad := AT3Zmad_freq / WT3Zmad_freq]

ggplot(residue_freq, aes(residue, log2(div10Z))) + geom_col()
ggplot(residue_freq, aes(residue, log2(div10Zmad))) + geom_col()
ggplot(residue_freq, aes(residue, log2(div3Z))) + geom_col()
ggplot(residue_freq, aes(residue, log2(div3Zmad))) + geom_col()

# Z-score using MAD correlates better to pause score
cor(ATzmad3A$AT1_A_pause, ATzmad3A$AT1_A_Zmad)
#[1] 0.8895451
cor(ATzmad3A$AT2_A_pause, ATzmad3A$AT2_A_Zmad)
#[1] 0.8816505
cor(ATzmad10A$AT1_A_pause, ATzmad10A$AT1_A_Zmad)
#[1] 0.8468532
cor(ATzmad10A$AT2_A_pause, ATzmad10A$AT2_A_Zmad)
#[1] 0.8356727
cor(ATz3A$AT1_A_pause, ATz3A$AT1_A_Z)
#[1] 0.6287343
cor(ATz3A$AT2_A_pause, ATz3A$AT2_A_Z)
#[1] 0.6261175
cor(ATz10A$AT1_A_pause, ATz10A$AT1_A_Z)
#[1] 0.5342286
cor(ATz10A$AT2_A_pause, ATz10A$AT2_A_Z)
#[1] 0.530547
cor(AT_P10$AT1_A_pause, ATzmad10A$AT1_A_Zmad)

### Calculate p-value using Fisher's exact test
his_fishers <- his_dt[WT1_A_rpc >= 0.5 & WT2_A_rpc >= 0.5 & 
                        AT1_A_rpc >= 0.5 & AT2_A_rpc >= 0.5 &
                        WT1_A_sum >= 64 & WT2_A_sum >= 64 &
                        AT1_A_sum >= 64 & AT2_A_sum >= 64]
View(his_fishers[(AT1_A_sum - AT1_A) < 0]) # Find orfs with position with aberrantly high number of reads
View(his_fishers[(AT2_A_sum - AT2_A) < 0])
View(his_fishers[(WT1_A_sum - WT1_A) < 0])
View(his_fishers[(WT2_A_sum - WT2_A) < 0])
his_fishers <- his_fishers[orf == "YCR012W"]
for (i in 1:nrow(his_fishers)) {
  counts1 <- matrix(c(his_fishers[i]$WT1_A, his_fishers[i]$AT1_A, 
                      (his_fishers[i]$WT1_A_sum - his_fishers[i]$WT1_A), 
                      (his_fishers[i]$AT1_A_sum - his_fishers[i]$AT1_A)), ncol = 2)
  counts2 <- matrix(c(his_fishers[i]$WT2_A, his_fishers[i]$AT2_A, 
                      (his_fishers[i]$WT2_A_sum - his_fishers[i]$WT2_A), 
                      (his_fishers[i]$AT2_A_sum - his_fishers[i]$AT2_A)), ncol = 2)
  stalling_test1 <- fisher.test(counts1)
  stalling_test2 <- fisher.test(counts2)
  his_fishers[i, pvalue1 := stalling_test1$p.value]
  his_fishers[i, pvalue2 := stalling_test2$p.value]
}
his_fishers[, padj1 := p.adjust(pvalue1, method = "BH"), by = orf]
his_fishers[, padj2 := p.adjust(pvalue2, method = "BH"), by = orf]
stalling_peaks <- his_fishers[AT1_A_Z > WT1_A_Z & AT2_A_Z > WT2_A_Z & AT1_A_Z > 3.5 & AT2_A_Z > 3.5 &
                                  padj1 < 0.05 & padj2 < 0.05]
View(stalling_peaks)
View(his_fishers[AT1_A_Z > WT1_A_Z & AT2_A_Z > WT2_A_Z & AT1_A_Z > 3.5 & AT2_A_Z > 3.5])

### Thoughts
# Can't use DESeq2 on genome basis because doesn't control for gene-specific differences, 
# but can't use it on gene-by-gene basis since not enough codons in smaller genes to be accurate.
# Fisher's might work depending on pvalue calculation for 2x3 table. However, seems unnecessary based statistical
# power provided by Z score (e.g. for PGK1, ommitting the pvalue threshold doesn't change the 
# identified positions if just use z-score).