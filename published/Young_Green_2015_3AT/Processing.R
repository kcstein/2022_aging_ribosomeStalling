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


