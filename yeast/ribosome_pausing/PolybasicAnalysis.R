library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
sc_dtA <- readRDS("sc_dtA.rds")

### Align occupancy at start of polybasic region ###
dir <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Yeast/PolybasicRegions/"
files <- list.files(dir, pattern = "*_only.csv")
samples <- sub("\\_only.csv", '', files)

for(i in 1:length(files)) {
  print(files[i])
  j <- as.data.table(lapply(paste0(dir, files[i]), read.csv, header = TRUE, stringsAsFactors = TRUE))
  setnames(j, c("orf", "polybasic", "name"))
  setkeyv(j, c("orf"))
  assign(samples[i], j)
}

KR3of3_dt <- sc_dtA[KR3of3, allow.cartesian = TRUE]
KR4of4_dt <- sc_dtA[KR4of4, allow.cartesian = TRUE]
KR5of5_dt <- sc_dtA[KR5of5, allow.cartesian = TRUE]
KR6of6_dt <- sc_dtA[KR6of6, allow.cartesian = TRUE]
KR3of3_dt[, adjusted := position - polybasic]
KR4of4_dt[, adjusted := position - polybasic]
KR5of5_dt[, adjusted := position - polybasic]
KR6of6_dt[, adjusted := position - polybasic]


### Normalize across polybasic interval ###
test <- KR6of6_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("name"))
setkeyv(KR6of6_dt, c("name"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = name]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = name]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = name]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = name]
test[, sch9_D0_1A_rpc_adjusted := mean(sch9_D0_1A), by = name]
test[, sch9_D0_2A_rpc_adjusted := mean(sch9_D0_2A), by = name]
test[, sch9_D4_1A_rpc_adjusted := mean(sch9_D4_1A), by = name]
test[, sch9_D4_2A_rpc_adjusted := mean(sch9_D4_2A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(KR6of6_dt$name, test$name))
KR6of6_dt <- cbind(KR6of6_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, sch9_D0_1A_rpc_adjusted = test[i]$sch9_D0_1A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, sch9_D0_2A_rpc_adjusted = test[i]$sch9_D0_2A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, sch9_D4_1A_rpc_adjusted = test[i]$sch9_D4_1A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, sch9_D4_2A_rpc_adjusted = test[i]$sch9_D4_2A_rpc_adjusted)
KR6of6_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
KR6of6_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
KR6of6_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
KR6of6_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
KR6of6_dt[, sch9_D0_1A_norm := sch9_D0_1A / sch9_D0_1A_rpc_adjusted]
KR6of6_dt[, sch9_D0_2A_norm := sch9_D0_2A / sch9_D0_2A_rpc_adjusted]
KR6of6_dt[, sch9_D4_1A_norm := sch9_D4_1A / sch9_D4_1A_rpc_adjusted]
KR6of6_dt[, sch9_D4_2A_norm := sch9_D4_2A / sch9_D4_2A_rpc_adjusted]
# saveRDS(KR6of6_dt, "KR6of6_dt.rds")


### Plot ribosome occupancy around region of interest ###
KR3of3_dt <- readRDS("doc/KR3of3_dt.rds")
KR4of4_dt <- readRDS("doc/KR4of4_dt.rds")
KR5of5_dt <- readRDS("doc/KR5of5_dt.rds")
KR6of6_dt <- readRDS("doc/KR6of6_dt.rds")

# average replicates and adjusted rpc0.1
ggplot(data = KR4of10_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & position > 25,
                         .(adjusted, D0A_pause = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=5, center=T),
                           D4A_pause = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=5, center=T))]) + xlim(-40, 40) +
  #stat_summary(aes(adjusted, D1normA), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
  #            fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(adjusted, D0A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  #stat_summary(aes(adjusted, D12normA), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
  #            fun.args=list(conf.int=0.5), fill = 'red') +
  stat_summary(aes(adjusted, D4A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'red')

# Sch9
ggplot(data = KR6of6_dt[sch9_D0_1A_rpc_adjusted >= 0.1 & sch9_D0_2A_rpc_adjusted >= 0.1 & 
                          sch9_D2_1A_rpc_adjusted >= 0.1 & sch9_D2_2A_rpc_adjusted >= 0.1 &
                          sch9_D4_1A_rpc_adjusted >= 0.1 & sch9_D4_2A_rpc_adjusted >= 0.1 & position > 25,
                        .(adjusted, D0_1A_pause = movingAverage(sch9_D0_1A_norm, n=5, center=T),
                          D0_2A_pause = movingAverage(sch9_D0_2A_norm, n=5, center=T),
                          D2_1A_pause = movingAverage(sch9_D2_1A_norm, n=5, center=T),
                          D2_2A_pause = movingAverage(sch9_D2_2A_norm, n=5, center=T),
                          D4_1A_pause = movingAverage(sch9_D4_1A_norm, n=5, center=T),
                          D4_2A_pause = movingAverage(sch9_D4_2A_norm, n=5, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D0_1A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D0_2A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D2_1A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(aes(adjusted, D2_2A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(aes(adjusted, D4_1A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(adjusted, D4_2A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'red')

ggplot(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                          WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 &
                          sch9_D0_1A_rpc_adjusted >= 0.1 & sch9_D0_2A_rpc_adjusted >= 0.1 & 
                          sch9_D4_1A_rpc_adjusted >= 0.1 & sch9_D4_2A_rpc_adjusted >= 0.1 & polybasic > 25 & stopdist < -25,
                        .(adjusted, D0A_pause = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=5, center=T),
                          D4A_pause = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=5, center=T),
                          S0A_pause = movingAverage(((sch9_D0_1A_norm + sch9_D0_2A_norm) / 2), n=5, center=T),
                          S4A_pause = movingAverage(((sch9_D4_1A_norm + sch9_D4_2A_norm) / 2), n=5, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D0A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D4A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(adjusted, S0A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(aes(adjusted, S4A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'orange')


ggplot(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                          WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 &
                          WT_D2_1A_rpc_adjusted >= 0.1 & WT_D2_2A_rpc_adjusted >= 0.1 & polybasic > 25 & stopdist < -25,
                        .(adjusted, D0A_pause = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                          D4A_pause = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T),
                          D2A_pause = movingAverage(((WT_D2_1A_norm + WT_D2_2A_norm) / 2), n=1, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D0A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D4A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(adjusted, D2A_pause), fun.y = "mean", geom = "line", size=0.5, color = 'blue')

ggplot(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                          WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 &
                          WT_D2_1A_rpc_adjusted >= 0.1 & WT_D2_2A_rpc_adjusted >= 0.1 & polybasic > 25 & stopdist < -25,
                        .(adjusted, D0_1 = movingAverage(WT_D0_1A_norm, n=5, center=T),
                          D0_2 = movingAverage(WT_D0_2A_norm, n=5, center=T),
                          D2_1 = movingAverage(WT_D2_1A_norm, n=5, center=T),
                          D2_2 = movingAverage(WT_D2_2A_norm, n=5, center=T),
                          D4_1 = movingAverage(WT_D4_1A_norm, n=5, center=T),
                          D4_2 = movingAverage(WT_D4_2A_norm, n=5, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D0_1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(adjusted, D0_2), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D2_1), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(aes(adjusted, D2_2), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(aes(adjusted, D4_1), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(adjusted, D4_2), fun.y = "mean", geom = "line", size=0.5, color = 'red')
