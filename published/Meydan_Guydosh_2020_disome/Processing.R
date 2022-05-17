library(data.table)
library(dplyr)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')

### Import codon density files ###
dir <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Meydan_Guydosh_2020/"
files <- list.files(dir, pattern = "*.codon")
samples <- sub("\\.codon", '', files)

for(i in 1:length(files)) {
  print(files[i])
  assign(samples[i], as.data.frame(lapply(paste0(dir, files[i]), read.delim, header = FALSE, stringsAsFactors = TRUE)))
}
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))


### Create data table with basic ORF properties ###
disome_dt <- data.table(orf = dfs[[1]][, 1],
                    position = dfs[[1]][, 2],
                    codon = dfs[[1]][, 3]) 
setkeyv(disome_dt, c("orf"))
orfs <- unique(disome_dt[, 1])
disome_dt[, length := (length(position) - 15), by = orf] # subtract flanking 7 codons and stop codon
disome_dt[, stopdist := position - (length + 1)] # stop codon is at stopdist == 0
disome_dt[, ID := as.character(paste(orf, position, sep = "_"))]

# Add residue
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(disome_dt$codon, codon_table$codon))
disome_dt <- cbind(disome_dt, residue = codon_table[i]$residue)

# Add raw counts
for(i in 1:length(samples)) {
  print(samples[i])
  disome_dt <- cbind(disome_dt, j = dfs[[i]][,4])
  setnames(disome_dt, "j", samples[i])
}


### Calculate rpc for each orf after excluding first and last 30 sense codons (Rachel Green excludes 100nt in the eIF5A paper)
# Sum reads in order to exclude genes with less than 64 reads in each replicate in downstream analysis
temp <- disome_dt[position > 20 & stopdist < -20]
for(i in samples) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  gene <- cbind(match(disome_dt$orf, temp$orf))
  disome_dt <- cbind(disome_dt, sum1 = temp[gene]$sum1)
  disome_dt <- cbind(disome_dt, rpc1 = temp[gene]$rpc1)
  setnames(disome_dt, "rpc1", paste0(i,"_rpc"))
  setnames(disome_dt, "sum1", paste0(i,"_sum"))
}
saveRDS(disome_dt, "disome_dt.rds")


### Calculate rpm, pause score ###
for(i in samples) {
  print(i)
  disome_dt[, paste0(i,"_rpm") := (disome_dt[[i]] / sum(disome_dt[[i]])) * 10^6] # calculate rpm
  disome_dt[, paste0(i,"_pause") := disome_dt[[i]] / disome_dt[[paste0(i,"_rpc")]]] # calculate pause
}
#saveRDS(disome_dt, "disome_dt.rds")

### Add 3 residue motif in active site ###
sc_dtA <- readRDS("../../Chronological/doc/sc_dtA.rds")
i <- cbind(match(disome_dt$ID, sc_dtA$ID))
disome_dt <- cbind(disome_dt, motif3 = sc_dtA[i]$motif3)
#saveRDS(disome_dt, "disome_dt.rds")

disome_dt[, WT_mono := (WT_mono1_A_pause + WT_mono2_A_pause) / 2]
disome_dt[, WT_di := (WT_di1_A_pause + WT_di2_A_pause + WT_di3_A_pause) / 3]
#saveRDS(disome_dt, "disome_dt.rds")
