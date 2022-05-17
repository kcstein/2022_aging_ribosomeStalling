library(data.table)
library(dplyr)


### Import codon density files ###
dir1 <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCSJF05_Ce/"
files1 <- list.files(dir1, pattern = "*.codon")
samples1 <- sub("\\.codon", '', files1)

for(i in 1:length(files1)) {
  print(files1[i])
  assign(samples1[i], as.data.frame(lapply(paste0(dir1, files1[i]), read.delim, header = FALSE, stringsAsFactors = TRUE)))
}

dir2a <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCSJF08_Ce_daf/"
files2a <- list.files(dir2a, pattern = "*.codon")
samples2a <- sub("\\.codon", '', files2a)
for(i in 1:length(files2a)) {
  print(files2a[i])
  assign(samples2a[i], as.data.frame(lapply(paste0(dir2a, files2a[i]), read.delim, header = FALSE, stringsAsFactors = TRUE)))
}

dir2b <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCSJF09_Ce_daf/"
files2b <- list.files(dir2b, pattern = "*.codon")
samples2b <- sub("\\.codon", '', files2b)
for(i in 1:length(files2b)) {
  print(files2b[i])
  assign(samples2b[i], as.data.frame(lapply(paste0(dir2b, files2b[i]), read.delim, header = FALSE, stringsAsFactors = TRUE)))
}

test <- append(files1, files2a)
files <- append(test, files2b)
samples <- sub("\\.codon", '', files)

dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))


### Create data table with basic ORF properties ###
ce_dt <- data.table(orf = dfs[[1]][, 1],
                    position = dfs[[1]][, 2],
                    codon = dfs[[1]][, 3]) 
setkeyv(ce_dt, c("orf"))
orfs <- unique(ce_dt[, 1])
ce_dt[, length := (length(position) - 15), by = orf] # subtract flanking 7 codons and stop codon
ce_dt[, stopdist := position - (length + 1)] # stop codon is at stopdist == 0
ce_dt[, ID := as.character(paste(orf, position, sep = "_"))]

# Add residue
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(ce_dt$codon, codon_table$codon))
ce_dt <- cbind(ce_dt, residue = codon_table[i]$residue)

# Add raw counts
for(i in 1:length(samples)) {
  print(samples[i])
  ce_dt <- cbind(ce_dt, j = dfs[[i]][,4])
  setnames(ce_dt, "j", samples[i])
}
#saveRDS(ce_dt, "ce_dt.rds")

ce_dt <- readRDS("ce_dt.rds")
N2samples <- samples[c(1:12,31:39,58:66)]
N2samplesA <- samples[c(1,4,7,10,58,61,64)]
N2_dt <- ce_dt[, c(1:19,38:46,65:73)]
N2_dtA <- ce_dt[, c(1:8,11,14,17,65,68,71)]
#saveRDS(N2_dt, "N2_dt.rds")
#saveRDS(N2_dtA, "N2_dtA.rds")
daf_samples <- samples[c(13:30,40:57)]
daf_dt <- ce_dt[, c(1:7,20:37,47:64)]
#saveRDS(daf_dt, "daf_dt.rds")


### Subset to just A-site alignment ###
N2_dtA <- N2_dt[, c(1:8,11,14,17,20,23,26,29,32,35)]
#N2samplesA <- names(N2_dtA[, c(8:17)])
N2samplesA <- names(N2_dtA[, c(8:11)])
saveRDS(N2_dtA, "N2_dtA.rds")

daf_dtA <- ce_dt[, c(1:7,20,23,26,29,32,35,47,50,53,56,59,62)]
daf_samplesA <- names(daf_dtA[, c(8:19)])
saveRDS(daf_dtA, "daf_dtA.rds")


### Calculate rpc for each orf after excluding first and last 30 sense codons (Rachel Green excludes 100nt in the eIF5A paper)
# Sum reads in order to exclude genes with less than 64 reads in each replicate in downstream analysis
N2_dt <- readRDS("N2_dt.rds")
temp <- N2_dt[position > 30 & stopdist < -30]
for(i in N2samples) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  gene <- cbind(match(N2_dt$orf, temp$orf))
  N2_dt <- cbind(N2_dt, sum1 = temp[gene]$sum1)
  N2_dt <- cbind(N2_dt, rpc1 = temp[gene]$rpc1)
  setnames(N2_dt, "rpc1", paste0(i,"_rpc"))
  setnames(N2_dt, "sum1", paste0(i,"_sum"))
}
saveRDS(N2_dt, "N2_dt.rds")

daf_dt <- readRDS("daf_dt.rds")
temp <- daf_dt[position > 30 & stopdist < -30]
for(i in daf_samples) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  temp[, med1 := lapply(.SD[, ..i], median), by = orf]
  gene <- cbind(match(daf_dt$orf, temp$orf))
  daf_dt <- cbind(daf_dt, sum1 = temp[gene]$sum1)
  daf_dt <- cbind(daf_dt, rpc1 = temp[gene]$rpc1)
  daf_dt <- cbind(daf_dt, med1 = temp[gene]$med1)
  setnames(daf_dt, "rpc1", paste0(i,"_rpc"))
  setnames(daf_dt, "sum1", paste0(i,"_sum"))
}
saveRDS(daf_dt, "daf_dt.rds")

N2_dtA <- N2_dtA[, c(1:14)]
N2samplesA <- names(N2_dtA[, c(8:14)])
N2_dtA <- readRDS("N2_dtA.rds")
temp <- N2_dtA[position > 20 & stopdist < -20]
for(i in N2samplesA) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, coverage1 := length(which(.SD[, ..i] != 0)) / (length - 40), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  temp[, sd1 := lapply(.SD[, ..i], sd), by = orf]
  temp[, Z1 := scale(.SD[, ..i], center = TRUE, scale = TRUE), by = orf]
  #temp[, med1 := lapply(.SD[, ..i], median), by = orf]
  gene <- cbind(match(N2_dtA$orf, temp$orf))
  N2_dtA <- cbind(N2_dtA, sum1 = temp[gene]$sum1)
  N2_dtA <- cbind(N2_dtA, coverage1 = temp[gene]$coverage1)
  N2_dtA <- cbind(N2_dtA, rpc1 = temp[gene]$rpc1)
  N2_dtA <- cbind(N2_dtA, sd1 = temp[gene]$sd1)
  ID1 <- cbind(match(N2_dtA$ID, temp$ID))
  N2_dtA <- cbind(N2_dtA, Z1 = temp[ID1]$Z1)
  #N2_dtA <- cbind(N2_dtA, med1 = temp[gene]$med1)
  setnames(N2_dtA, "sum1", paste0(i,"_sum"))
  setnames(N2_dtA, "coverage1", paste0(i,"_coverage"))
  setnames(N2_dtA, "rpc1", paste0(i,"_rpc"))
  setnames(N2_dtA, "sd1", paste0(i,"_sd"))
  setnames(N2_dtA, "Z1", paste0(i,"_Z"))
  #setnames(N2_dtA, "med1", paste0(i,"_med"))
}
saveRDS(N2_dtA, "N2_dtA.rds")

daf_dtA <- readRDS("daf_dtA.rds")
temp <- daf_dtA[position > 30 & stopdist < -30]
for(i in daf_samplesA) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  temp[, med1 := lapply(.SD[, ..i], median), by = orf]
  gene <- cbind(match(daf_dtA$orf, temp$orf))
  daf_dtA <- cbind(daf_dtA, sum1 = temp[gene]$sum1)
  daf_dtA <- cbind(daf_dtA, rpc1 = temp[gene]$rpc1)
  daf_dtA <- cbind(daf_dtA, med1 = temp[gene]$med1)
  setnames(daf_dtA, "rpc1", paste0(i,"_rpc"))
  setnames(daf_dtA, "sum1", paste0(i,"_sum"))
  setnames(daf_dtA, "med1", paste0(i,"_med"))
}
saveRDS(daf_dtA, "daf_dtA.rds")


### Calculate rpm, pause score, and z-score ###
N2_dt <- readRDS("N2_dt.rds")
for(i in N2samples) {
  print(i)
  N2_dt[, paste0(i,"_rpm") := (N2_dt[[i]] / sum(N2_dt[[i]])) * 10^6] # calculate rpm
  N2_dt[, paste0(i,"_pause") := N2_dt[[i]] / N2_dt[[paste0(i,"_rpc")]]] # calculate pause
  N2_dt[, Z := abs(N2_dt[[i]] - N2_dt[[paste0(i,"_rpc")]])]
  N2_dt[, Z := median(Z), by = orf]
  N2_dt[, Z  := (0.6745*(N2_dt[[i]] - N2_dt[[paste0(i,"_rpc")]])) / Z] # calculate z-score
  setnames(N2_dt, "Z", paste0(i,"_Z"))
}
saveRDS(N2_dt, "N2_dt.rds")

daf_dt <- readRDS("daf_dt.rds")
for(i in daf_samples) {
  print(i)
  daf_dt[, paste0(i,"_rpm") := (daf_dt[[i]] / sum(daf_dt[[i]])) * 10^6] # calculate rpm
  daf_dt[, paste0(i,"_pause") := daf_dt[[i]] / daf_dt[[paste0(i,"_rpc")]]] # calculate pause
  daf_dt[, Z := abs(daf_dt[[i]] - daf_dt[[paste0(i,"_rpc")]])]
  daf_dt[, Z := median(Z), by = orf]
  daf_dt[, Z  := (0.6745*(daf_dt[[i]] - daf_dt[[paste0(i,"_rpc")]])) / Z] # calculate z-score
  setnames(daf_dt, "Z", paste0(i,"_Z"))
}
saveRDS(daf_dt, "daf_dt.rds")

for(i in N2samplesA) {
  print(i)
  N2_dtA[, paste0(i,"_rpm") := (N2_dtA[[i]] / sum(N2_dtA[[i]])) * 10^6] # calculate rpm
  N2_dtA[, paste0(i,"_pause") := N2_dtA[[i]] / N2_dtA[[paste0(i,"_rpc")]]] # calculate pause
  #N2_dtA[, Z := abs(N2_dtA[[i]] - N2_dtA[[paste0(i,"_med")]])]
  #N2_dtA[, Z := median(Z), by = orf]
  #N2_dtA[, Z := (0.6745*(N2_dtA[[i]] - N2_dtA[[paste0(i,"_med")]])) / Z] # calculate z-score
  #setnames(N2_dtA, "Z", paste0(i,"_Z"))
  #N2_dtA[, paste0(i,"_Z") := (N2_dtA[[i]] - N2_dtA[[paste0(i,"_rpc")]]) / ((N2_dtA[[paste0(i,"_rpc")]])^(1/2))]
}
saveRDS(N2_dtA, "N2_dtA.rds")

for(i in daf_samplesA) {
  print(i)
  daf_dtA[, paste0(i,"_rpm") := (daf_dtA[[i]] / sum(daf_dtA[[i]])) * 10^6] # calculate rpm
  daf_dtA[, paste0(i,"_pause") := daf_dtA[[i]] / daf_dtA[[paste0(i,"_rpc")]]] # calculate pause
  daf_dtA[, Z := abs(daf_dtA[[i]] - daf_dtA[[paste0(i,"_med")]])]
  daf_dtA[, Z := median(Z), by = orf]
  daf_dtA[, Z := (0.6745*(daf_dtA[[i]] - daf_dtA[[paste0(i,"_med")]])) / Z] # calculate z-score
  setnames(daf_dtA, "Z", paste0(i,"_Zmad"))
  daf_dtA[, paste0(i,"_Z") := (daf_dtA[[i]] - daf_dtA[[paste0(i,"_rpc")]]) / ((daf_dtA[[paste0(i,"_rpc")]])^(1/2))]
}
saveRDS(daf_dtA, "daf_dtA.rds")


### Add motif in active site ###
temp1 <- N2_dtA[position >= -1 & stopdist <= -1, c(1,2,5:7)]
#orfs <- unique(temp1[, 1])
#match(c("ZK973.8_1"), name)
#name1 <- name[1:496507]
#motif1 <- motif[1:496507]
dt <- data.table(ID = name1, motif3 = motif1)
saveRDS(dt, "dt9.rds")
#which(orfs$orf %in% c("ZK973.8"))
temp1 <- readRDS("temp1.rds")
orfs <- orfs[c(1150:1157)]
temp1 <- temp1[temp1$orf %in% orfs$orf]
#saveRDS(temp1, "temp1.rds")
name <- NULL
motif <- NULL
for (i in orfs$orf) {
  print(i)
  temp <- temp1[orf == i]
  for (j in 3:nrow(temp)) {
    name1 <- temp[j]$ID
    name <- c(name, name1)
    E <- temp[j-2]$residue
    P <- temp[j-1]$residue
    A <- temp[j]$residue
    motif3 <- paste0(E,P,A)
    motif <- c(motif, motif3)
  }
}
dt <- data.table(ID = name, motif3 = motif)
dt1 <- readRDS("dt1.rds")
dt2 <- readRDS("dt2.rds")
dt3 <- readRDS("dt3.rds")
dt4 <- readRDS("dt4.rds")
dt5 <- readRDS("dt5.rds")
dt6 <- readRDS("dt6.rds")
dt7 <- readRDS("dt7.rds")
dt8 <- readRDS("dt8.rds")
dt9 <- readRDS("dt9.rds")
new_dt <- rbind(dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt)
i <- cbind(match(N2_dtA$ID, new_dt$ID))
N2_dtA <- cbind(N2_dtA, motif3 = new_dt[i]$motif3)
N2_dtA[, motif2 := substring(motif3,2,3)]
N2_dtA[, D1_pause := (D1_1A_pause + D1_2A_pause) / 2]
N2_dtA[, D12_pause := (D12_1A_pause + D12_2A_pause) / 2]
saveRDS(N2_dtA, "N2_dtA.rds")

names(N2_dtA)
N2_dtA <- N2_dtA[, c(1:11,78,81)]


N2_dtA_archive <- readRDS("doc/Archive_Misc/N2_dtA_2019.09.05.rds")
i <- cbind(match(N2_dtA$ID, N2_dtA_archive$ID))
N2_dtA <- cbind(N2_dtA, motif3 = N2_dtA_archive[i]$motif3)
N2_dtA[, motif2 := substring(motif3,2,3)]
N2_dtA[, D1_pause := (D1_1A_pause + D1_2A_pause) / 2]
N2_dtA[, D12_pause := (D12_1A_pause + D12_2A_pause) / 2]
N2_dtA[, D1_pause_all := (D1_1A_pause + D1_2A_pause + N2_D1_2A_pause) / 3]
N2_dtA[, D12_pause_all := (D12_1A_pause + D12_2A_pause + N2_D12_2A_pause) / 3]
saveRDS(N2_dtA, "N2_dtA.rds")


### Add gene names
genenames <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/geneNames.csv', header = TRUE, stringsAsFactors = TRUE))
hartl_proteome <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/proteinLists/hartl_proteome.csv', header = TRUE, stringsAsFactors = TRUE))
i <- cbind(match(N2_dtA$orf, genenames$Transcript.stable.ID))
N2_dtA <- cbind(N2_dtA, WBgene = genenames[i]$Gene.stable.ID)
i <- cbind(match(N2_dtA$WBgene, hartl_proteome$WBgene))
N2_dtA <- cbind(N2_dtA, gene = hartl_proteome[i]$gene)
N2_dtA <- cbind(N2_dtA, uniprot = hartl_proteome[i]$uniprot)
orfs <- N2_dtA[, .SD[which.max(position)], by = orf]
saveRDS(N2_dtA, "N2_dtA.rds")

