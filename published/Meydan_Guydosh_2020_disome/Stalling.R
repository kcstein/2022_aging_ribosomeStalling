library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
disome_dt <- readRDS("disome_dt.rds")

### Identify max disome position in each gene ###
disome_positions <- disome_dt[, .SD[which.max(WT_di)], by = orf]
disome_positions[, ratio := WT_di / WT_mono]
disome_positions <- disome_positions[ratio > 2 & WT_di < Inf & ratio < Inf]
write.csv(disome_positions[, c(1,2,6,46)], "/Users/KevinStein/Desktop/disome_positions.csv")
