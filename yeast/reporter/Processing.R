library(data.table)
library(dplyr)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')

### Import codon density files ###
dir <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCSJF06_reporters/"
files <- list.files(dir, pattern = "*.codon")
samples <- sub("\\.codon", '', files)

for(i in 1:length(files)) {
  print(files[i])
  assign(samples[i], as.data.frame(lapply(paste0(dir, files[i]), read.delim, header = FALSE, stringsAsFactors = TRUE)))
}

### Create data table with basic ORF properties ###
GH_dt <- data.table(orf = GH0_1A[, 1],
                        position = GH0_1A[, 2],
                        codon = GH0_1A[, 3]) 
setkeyv(GH_dt, c("orf"))
orfs <- unique(GH_dt[, 1])
GH_dt[, length := (length(position) - 15), by = orf] # subtract flanking 7 codons and stop codon
GH_dt[, stopdist := position - (length + 1)] # stop codon is at stopdist == 0
GH_dt[, ID := as.character(paste(orf, position, sep = "_"))]

# Add residue
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(GH_dt$codon, codon_table$codon))
GH_dt <- cbind(GH_dt, residue = codon_table[i]$residue)

# Add raw counts
GH_dt <- cbind(GH_dt, GH0_1A = GH0_1A[, 4])
GH_dt <- cbind(GH_dt, GH0_2A = GH0_2A[, 4])
GH_dt <- cbind(GH_dt, GH4_1A = GH4_1A[, 4])
GH_dt <- cbind(GH_dt, GH4_2A = GH4_2A[, 4])


### Calculate rpc for each orf after excluding first and last 30 sense codons (Rachel Green excludes 100nt in the eIF5A paper)
# Sum reads in order to exclude genes with less than 64 reads in each replicate in downstream analysis
temp <- GH_dt[position > 20 & stopdist < -20]
GHsamples <- samples[1:4]
for(i in GHsamples) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  gene <- cbind(match(GH_dt$orf, temp$orf))
  GH_dt <- cbind(GH_dt, sum1 = temp[gene]$sum1)
  GH_dt <- cbind(GH_dt, rpc1 = temp[gene]$rpc1)
  setnames(GH_dt, "rpc1", paste0(i,"_rpc"))
  setnames(GH_dt, "sum1", paste0(i,"_sum"))
}

### Calculate rpm, pause score ###
for(i in GHsamples) {
  print(i)
  GH_dt[, paste0(i,"_rpm") := (GH_dt[[i]] / sum(GH_dt[[i]])) * 10^6] # calculate rpm
  GH_dt[, paste0(i,"_pause") := GH_dt[[i]] / GH_dt[[paste0(i,"_rpc")]]] # calculate pause
}
#saveRDS(GH_dt, "GH_dt.rds")

### Add 3 residue motif in active site ###
sc_dtA <- readRDS("../Chronological/doc/sc_dtA.rds")
i <- cbind(match(GH_dt$ID, sc_dtA$ID))
GH_dt <- cbind(GH_dt, motif3 = sc_dtA[i]$motif3)

GH_dt[, GH0_pause := (GH0_1A_pause + GH0_2A_pause) / 2]
GH_dt[, GH4_pause := (GH4_1A_pause + GH4_2A_pause) / 2]
#saveRDS(GH_dt, "GH_dt.rds")


### Create data table with basic ORF properties ###
GHA_dt <- data.table(orf = GHA0_1A[, 1],
                     position = GHA0_1A[, 2],
                     codon = GHA0_1A[, 3]) 
setkeyv(GHA_dt, c("orf"))
orfs <- unique(GHA_dt[, 1])
GHA_dt[, length := (length(position) - 15), by = orf] # subtract flanking 7 codons and stop codon
GHA_dt[, stopdist := position - (length + 1)] # stop codon is at stopdist == 0
GHA_dt[, ID := as.character(paste(orf, position, sep = "_"))]

# Add residue
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(GHA_dt$codon, codon_table$codon))
GHA_dt <- cbind(GHA_dt, residue = codon_table[i]$residue)

# Add raw counts
GHA_dt <- cbind(GHA_dt, GHA0_1A = GHA0_1A[, 4])
GHA_dt <- cbind(GHA_dt, GHA0_2A = GHA0_2A[, 4])
GHA_dt <- cbind(GHA_dt, GHA4_1A = GHA4_1A[, 4])
GHA_dt <- cbind(GHA_dt, GHA4_2A = GHA4_2A[, 4])


### Calculate rpc for each orf after excluding first and last 30 sense codons (Rachel Green excludes 100nt in the eIF5A paper)
# Sum reads in order to exclude genes with less than 64 reads in each replicate in downstream analysis
temp <- GHA_dt[position > 20 & stopdist < -20]
GHAsamples <- samples[5:8]
for(i in GHAsamples) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  gene <- cbind(match(GHA_dt$orf, temp$orf))
  GHA_dt <- cbind(GHA_dt, sum1 = temp[gene]$sum1)
  GHA_dt <- cbind(GHA_dt, rpc1 = temp[gene]$rpc1)
  setnames(GHA_dt, "rpc1", paste0(i,"_rpc"))
  setnames(GHA_dt, "sum1", paste0(i,"_sum"))
}

### Calculate rpm, pause score ###
for(i in GHAsamples) {
  print(i)
  GHA_dt[, paste0(i,"_rpm") := (GHA_dt[[i]] / sum(GHA_dt[[i]])) * 10^6] # calculate rpm
  GHA_dt[, paste0(i,"_pause") := GHA_dt[[i]] / GHA_dt[[paste0(i,"_rpc")]]] # calculate pause
}

### Add 3 residue motif in active site ###
sc_dtA <- readRDS("../Chronological/doc/sc_dtA.rds")
i <- cbind(match(GHA_dt$ID, sc_dtA$ID))
GHA_dt <- cbind(GHA_dt, motif3 = sc_dtA[i]$motif3)

GHA_dt[, GHA0_pause := (GHA0_1A_pause + GHA0_2A_pause) / 2]
GHA_dt[, GHA4_pause := (GHA4_1A_pause + GHA4_2A_pause) / 2]
#saveRDS(GHA_dt, "GHA_dt.rds")

