library(data.table)
library(dplyr)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')

### Import codon density files ###
dir <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCSJF10_CA/UTR/"
files <- list.files(dir, pattern = "*.codon")
samples <- sub("\\.codon", '', files)

for(i in 1:length(files)) {
  print(files[i])
  assign(samples[i], as.data.frame(lapply(paste0(dir, files[i]), read.delim, header = FALSE, stringsAsFactors = TRUE)))
}
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))


### Create data table with basic ORF properties ###
sc_dt <- data.table(orf = dfs[[1]][, 1],
                    position = dfs[[1]][, 2],
                    codon = dfs[[1]][, 3]) 
setkeyv(sc_dt, c("orf"))
orfs <- unique(sc_dt[, 1])
sc_dt[, length := (length(position) - 15), by = orf] # subtract flanking 7 codons and stop codon
sc_dt[, stopdist := position - (length + 1)] # stop codon is at stopdist == 0
sc_dt[, ID := as.character(paste(orf, position, sep = "_"))]

scCR_dt <- data.table(orf = dfsCR[[1]][, 1],
                    position = dfsCR[[1]][, 2],
                    codon = dfsCR[[1]][, 3]) 
setkeyv(scCR_dt, c("orf"))
orfsCR <- unique(scCR_dt[, 1])
scCR_dt[, length := (length(position) - 15), by = orf] # subtract flanking 7 codons and stop codon
scCR_dt[, stopdist := position - (length + 1)] # stop codon is at stopdist == 0
scCR_dt[, ID := as.character(paste(orf, position, sep = "_"))]

# Add residue
codon_table <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/codon_table.csv"))
i <- cbind(match(sc_dt$codon, codon_table$codon))
sc_dt <- cbind(sc_dt, residue = codon_table[i]$residue)
i <- cbind(match(scCR_dt$codon, codon_table$codon))
scCR_dt <- cbind(scCR_dt, residue = codon_table[i]$residue)

# Add raw counts
for(i in 1:length(samples)) {
  print(samples[i])
  sc_dt <- cbind(sc_dt, j = dfs[[i]][,4])
  setnames(sc_dt, "j", samples[i])
}

for(i in 1:length(samplesCR)) {
  print(samplesCR[i])
  scCR_dt <- cbind(scCR_dt, j = dfsCR[[i]][,4])
  setnames(scCR_dt, "j", samplesCR[i])
}


### Subset to A-site alignment ###
sc_dtA <- sc_dt[, c(1:8,11,14,17,20,23,26,29,32,35,38,41)]
samplesA <- names(sc_dtA[, c(8:19)])
#saveRDS(sc_dtA, "sc_dtA.rds")


### Calculate rpc for each orf after excluding first and last 30 sense codons (Rachel Green excludes 100nt in the eIF5A paper)
# Sum reads in order to exclude genes with less than 64 reads in each replicate in downstream analysis
temp <- sc_dt[position > 20 & stopdist < -20]
for(i in samples) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  gene <- cbind(match(sc_dt$orf, temp$orf))
  sc_dt <- cbind(sc_dt, sum1 = temp[gene]$sum1)
  sc_dt <- cbind(sc_dt, rpc1 = temp[gene]$rpc1)
  setnames(sc_dt, "rpc1", paste0(i,"_rpc"))
  setnames(sc_dt, "sum1", paste0(i,"_sum"))
}
saveRDS(sc_dt, "sc_dt.rds")


sc_dtA <- sc_dtA[, c(1:21)]
temp <- sc_dtA[position > 20 & stopdist < -20]
for(i in samplesA) {
  print(i)
  temp[, sum1 := lapply(.SD[, ..i], sum), by = orf]
  temp[, coverage1 := length(which(.SD[, ..i] != 0)) / (length - 40), by = orf]
  temp[, rpc1 := lapply(.SD[, ..i], mean), by = orf]
  temp[, sd1 := lapply(.SD[, ..i], sd), by = orf]
  temp[, Z1 := scale(.SD[, ..i], center = TRUE, scale = TRUE), by = orf]
  #temp[, med1 := lapply(.SD[, ..i], median), by = orf]
  gene <- cbind(match(sc_dtA$orf, temp$orf))
  sc_dtA <- cbind(sc_dtA, sum1 = temp[gene]$sum1)
  sc_dtA <- cbind(sc_dtA, coverage1 = temp[gene]$coverage1)
  sc_dtA <- cbind(sc_dtA, rpc1 = temp[gene]$rpc1)
  sc_dtA <- cbind(sc_dtA, sd1 = temp[gene]$sd1)
  ID1 <- cbind(match(sc_dtA$ID, temp$ID))
  sc_dtA <- cbind(sc_dtA, Z1 = temp[ID1]$Z1)
  #sc_dtA <- cbind(sc_dtA, med1 = temp[gene]$med1)
  setnames(sc_dtA, "sum1", paste0(i,"_sum"))
  setnames(sc_dtA, "coverage1", paste0(i,"_coverage"))
  setnames(sc_dtA, "rpc1", paste0(i,"_rpc"))
  setnames(sc_dtA, "sd1", paste0(i,"_sd"))
  setnames(sc_dtA, "Z1", paste0(i,"_Z"))
  #setnames(sc_dtA, "med1", paste0(i,"_med"))
}
saveRDS(sc_dtA, "sc_dtA.rds")


### Calculate rpm, pause score, and z-score ###
sc_dt <- readRDS("sc_dt.rds")
for(i in samples) {
  print(i)
  sc_dt[, paste0(i,"_rpm") := (sc_dt[[i]] / sum(sc_dt[[i]])) * 10^6] # calculate rpm
  sc_dt[, paste0(i,"_pause") := sc_dt[[i]] / sc_dt[[paste0(i,"_rpc")]]] # calculate pause
  #sc_dt[, Z := abs(sc_dt[[i]] - sc_dt[[paste0(i,"_rpc")]])]
  #sc_dt[, Z := median(Z), by = orf]
  #sc_dt[, Z  := (0.6745*(sc_dt[[i]] - sc_dt[[paste0(i,"_rpc")]])) / Z] # calculate z-score
  #setnames(sc_dt, "Z", paste0(i,"_Z"))
}
saveRDS(sc_dt, "sc_dt.rds")


sc_dtA <- readRDS("sc_dtA.rds")
for(i in samplesA) {
  print(i)
  sc_dtA[, paste0(i,"_rpm") := (sc_dtA[[i]] / sum(sc_dtA[[i]])) * 10^6] # calculate rpm
  sc_dtA[, paste0(i,"_pause") := sc_dtA[[i]] / sc_dtA[[paste0(i,"_rpc")]]] # calculate pause
  #sc_dtA[, Z := abs(sc_dtA[[i]] - sc_dtA[[paste0(i,"_med")]])]
  #sc_dtA[, Z := median(Z), by = orf]
  #sc_dtA[, Z := (0.6745*(sc_dtA[[i]] - sc_dtA[[paste0(i,"_med")]])) / Z] # calculate z-score
  #setnames(sc_dtA, "Z", paste0(i,"_Z"))
}
saveRDS(sc_dtA, "sc_dtA.rds")


### Add 3 residue motif in active site ###
temp1 <- sc_dtA[position >= -1 & stopdist <= -1, c(1:3,6:7)]
#orfs <- unique(temp1[, 1])
#match(c("YOR347C_1"), name)
#name1 <- name[1:524057]
#motif1 <- motif[1:524057]
dt <- data.table(ID = name1, motif3 = motif1)
saveRDS(dt, "dt3.rds")
#which(orfs$orf %in% c("YOR347C"))
temp1 <- readRDS("temp1.rds")
orfs <- orfs[c(1090:1590)]
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
new_dt <- rbind(dt1,dt2,dt3,dt)
i <- cbind(match(sc_dtA$ID, new_dt$ID))
sc_dtA <- cbind(sc_dtA, motif3 = new_dt[i]$motif3)
sc_dtA[, motif2 := substring(motif3,2,3)]
sc_dtA[, WT_D0_pause := (WT_D0_1A_pause + WT_D0_2A_pause) / 2]
sc_dtA[, WT_D2_pause := (WT_D2_1A_pause + WT_D2_2A_pause) / 2]
sc_dtA[, WT_D4_pause := (WT_D4_1A_pause + WT_D4_2A_pause) / 2]
sc_dtA[, sch9_D0_pause := (sch9_D0_1A_pause + sch9_D0_2A_pause) / 2]
sc_dtA[, sch9_D2_pause := (sch9_D2_1A_pause + sch9_D2_2A_pause) / 2]
sc_dtA[, sch9_D4_pause := (sch9_D4_1A_pause + sch9_D4_2A_pause) / 2]
saveRDS(sc_dtA, "sc_dtA.rds")

names(sc_dtA)
sc_dtA <- sc_dtA[, c(1:19,92,97)]


### Add 2 codon motif in active site ###
temp1 <- sc_dtA[position >= -1 & stopdist <= -1, c(1:3,6)]
orfs <- unique(temp1[, 1])
#match(c("YOR347C_1"), name)
#name1 <- name[1:524057]
#motif1 <- motif[1:524057]
dt <- data.table(ID = name1, motif3 = motif1)
saveRDS(dt, "dt3.rds")
#which(orfs$orf %in% c("YOR347C"))
temp1 <- readRDS("temp1.rds")
orfs <- orfs[c(1090:1590)]
temp1 <- temp1[temp1$orf %in% orfs$orf]
#saveRDS(temp1, "temp1.rds")
name <- NULL
dicodon <- NULL
for (i in orfs$orf) {
  print(i)
  temp <- temp1[orf == i]
  for (j in 2:nrow(temp)) {
    name1 <- temp[j]$ID
    name <- c(name, name1)
    P <- temp[j-1]$codon
    A <- temp[j]$codon
    dicodon1 <- paste0(E,P,A)
    dicodon <- c(dicodon, dicodon1)
  }
}
dt <- data.table(ID = name, motif3 = motif)
dt1 <- readRDS("dt1.rds")
dt2 <- readRDS("dt2.rds")
dt3 <- readRDS("dt3.rds")
new_dt <- rbind(dt1,dt2,dt3,dt)
i <- cbind(match(sc_dtA$ID, new_dt$ID))
sc_dtA <- cbind(sc_dtA, motif3 = new_dt[i]$motif3)


CGA_dt <- sc_dtA[position >= 20 & stopdist <= -20 & codon == "CGA"]
CGA_dt[, CGA_1 := position - 1]
CGA_dt[, CGA_ID := as.character(paste(orf, CGA_1, sep = "_"))]
CGA_1_dt <- sc_dtA[sc_dtA$ID %in% CGA_dt$CGA_ID]
AGGCGA <- CGA_1_dt[codon == "AGG"]
ATACGA <- CGA_1_dt[codon == "ATA"]
CGACGA <- CGA_1_dt[codon == "CGA"]
CTGCGA <- CGA_1_dt[codon == "CTG"]
GTACGA <- CGA_1_dt[codon == "GTA"]
GTGCGA <- CGA_1_dt[codon == "GTG"]

CGG_dt <- sc_dtA[position >= 20 & stopdist <= -20 & codon == "CGG"]
CGG_dt[, CGG_1 := position - 1]
CGG_dt[, CGG_ID := as.character(paste(orf, CGG_1, sep = "_"))]
CGG_1_dt <- sc_dtA[sc_dtA$ID %in% CGG_dt$CGG_ID]
AGGCGG <- CGG_1_dt[codon == "AGG"]
ATACGG <- CGG_1_dt[codon == "ATA"]
CGACGG <- CGG_1_dt[codon == "CGA"]

CCG_dt <- sc_dtA[position >= 20 & stopdist <= -20 & codon == "CCG"]
CCG_dt[, CCG_1 := position - 1]
CCG_dt[, CCG_ID := as.character(paste(orf, CCG_1, sep = "_"))]
CCG_1_dt <- sc_dtA[sc_dtA$ID %in% CCG_dt$CCG_ID]
CGACCG <- CCG_1_dt[codon == "CGA"]
CTCCCG <- CCG_1_dt[codon == "CTC"]
CTGCCG <- CCG_1_dt[codon == "CTG"]
GTACCG <- CCG_1_dt[codon == "GTA"]

CTG_dt <- sc_dtA[position >= 20 & stopdist <= -20 & codon == "CTG"]
CTG_dt[, CTG_1 := position - 1]
CTG_dt[, CTG_ID := as.character(paste(orf, CTG_1, sep = "_"))]
CTG_1_dt <- sc_dtA[sc_dtA$ID %in% CTG_dt$CTG_ID]
CGACTG <- CTG_1_dt[codon == "CGA"]

GCG_dt <- sc_dtA[position >= 20 & stopdist <= -20 & codon == "GCG"]
GCG_dt[, GCG_1 := position - 1]
GCG_dt[, GCG_ID := as.character(paste(orf, GCG_1, sep = "_"))]
GCG_1_dt <- sc_dtA[sc_dtA$ID %in% CGG_dt$GCG_ID]
CGAGCG <- GCG_1_dt[codon == "CGA"]

ATA_dt <- sc_dtA[position >= 20 & stopdist <= -20 & codon == "ATA"]
ATA_dt[, ATA_1 := position - 1]
ATA_dt[, ATA_ID := as.character(paste(orf, ATA_1, sep = "_"))]
ATA_1_dt <- sc_dtA[sc_dtA$ID %in% ATA_dt$ATA_ID]
CGAATA <- ATA_1_dt[codon == "CGA"]
CTGATA <- ATA_1_dt[codon == "CTG"]

dicodon <- rbind(AGGCGA, AGGCGG, ATACGA, ATACGG, CGAATA, CGACCG, CGACGA, CGACGG, CGACTG, CGAGCG, CTCCCG, CTGATA, CTGCCG, CTGCGA, GTACCG, GTACGA, GTGCGA)
saveRDS(dicodon, "dicodon.rds")

dicodon_R <- rbind(AGGCGA, AGGCGG, ATACGA, ATACGG, CGACGA, CGACGG, CTGCGA, GTACGA, GTGCGA)


