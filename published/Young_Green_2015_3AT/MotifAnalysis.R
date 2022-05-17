library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
his_dt <- readRDS("doc/his_dt.rds")

### Polybasic analysis
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

his_dt_temp <- his_dt[, c(1:8,11,14,17,21,27,33,39)]
KR3of3_dt <- his_dt_temp[KR3of3, allow.cartesian=TRUE]
KR4of4_dt <- his_dt_temp[KR4of4, allow.cartesian=TRUE]
KR5of5_dt <- his_dt_temp[KR5of5, allow.cartesian=TRUE]
KR6of6_dt <- his_dt_temp[KR6of6, allow.cartesian=TRUE]
KR3of3_dt[, adjusted := position - polybasic]
KR4of4_dt[, adjusted := position - polybasic]
KR5of5_dt[, adjusted := position - polybasic]
KR6of6_dt[, adjusted := position - polybasic]


### Normalize across polybasic interval
test <- KR6of6_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("name"))
setkeyv(KR6of6_dt, c("name"))
test[, AT1_A_rpc_adjusted := mean(AT1_A), by = name]
test[, AT2_A_rpc_adjusted := mean(AT2_A), by = name]
test[, WT1_A_rpc_adjusted := mean(WT1_A), by = name]
test[, WT2_A_rpc_adjusted := mean(WT2_A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(KR6of6_dt$name, test$name))
KR6of6_dt <- cbind(KR6of6_dt, AT1_A_rpc_adjusted = test[i]$AT1_A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, AT2_A_rpc_adjusted = test[i]$AT2_A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT1_A_rpc_adjusted = test[i]$WT1_A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT2_A_rpc_adjusted = test[i]$WT2_A_rpc_adjusted)
KR6of6_dt[, AT1_A_norm := AT1_A / AT1_A_rpc_adjusted]
KR6of6_dt[, AT2_A_norm := AT2_A / AT2_A_rpc_adjusted]
KR6of6_dt[, WT1_A_norm := WT1_A / WT1_A_rpc_adjusted]
KR6of6_dt[, WT2_A_norm := WT2_A / WT2_A_rpc_adjusted]

# saveRDS(KR3of3_dt, "KR3of3_dt.rds")
# saveRDS(KR4of4_dt, "KR4of4_dt.rds")
# saveRDS(KR5of5_dt, "KR5of5_dt.rds")
# saveRDS(KR6of6_dt, "KR6of6_dt.rds")


### Plot ribosome occupancy around residue of interest
KR3of3_dt <- readRDS("doc/KR3of3_dt.rds")
KR4of4_dt <- readRDS("doc/KR4of4_dt.rds")
KR5of5_dt <- readRDS("doc/KR5of5_dt.rds")
KR6of6_dt <- readRDS("doc/KR6of6_dt.rds")

# average replicates, adjusted rpc0.1
plot <- ggplot(data = KR6of6_dt[WT1_A_rpc_adjusted >= 0.1 & WT2_A_rpc_adjusted >= 0.1,
                        .(adjusted, WT6 = movingAverage(((WT1_A_norm + WT2_A_norm) / 2), n=3, center=T))]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT6), fun.y = "mean", geom = "line", size=1.25, color = 'red') + 
  stat_summary(data = KR5of5_dt[WT1_A_rpc_adjusted >= 0.1 & WT2_A_rpc_adjusted >= 0.1,
            .(adjusted, WT5 = movingAverage(((WT1_A_norm + WT2_A_norm) / 2), n=3, center=T))], 
            aes(adjusted, WT5), fun.y = "mean", geom = "line", size=1.25, color = 'orange') + 
  stat_summary(data = KR4of4_dt[WT1_A_rpc_adjusted >= 0.1 & WT2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT4 = movingAverage(((WT1_A_norm + WT2_A_norm) / 2), n=3, center=T))], 
               aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'blue') + 
  stat_summary(data = KR3of3_dt[WT1_A_rpc_adjusted >= 0.1 & WT2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT3 = movingAverage(((WT1_A_norm + WT2_A_norm) / 2), n=3, center=T))], 
               aes(adjusted, WT3), fun.y = "mean", geom = "line", size=1.25, color = 'black')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Green_polybasic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### Examples
ggplot(data=his_dt[orf == "YAL010C", .(position, WT = movingAverage(((WT1_A_rpm + WT2_A_rpm) / 2), n=5, center=T), 
                                       AT = movingAverage(((AT1_A_rpm + AT2_A_rpm) / 2), n=5, center=T))]) + 
  geom_vline(xintercept = 44, color = "gray50", linetype = "dashed") +
  geom_vline(xintercept = 132, color = "gray50", linetype = "dashed") +
  geom_vline(xintercept = 277, color = "gray50", linetype = "dashed") +
  geom_vline(xintercept = 351, color = "gray50", linetype = "dashed") +
  geom_vline(xintercept = 374, color = "gray50", linetype = "dashed") +
  # geom_vline(xintercept = 209, color = "gray50", linetype = "dashed") +
  # geom_vline(xintercept = 248, color = "gray50", linetype = "dashed") +
  # geom_vline(xintercept = 300, color = "gray50", linetype = "dashed") +
  geom_line(aes(position, AT), color = "red") +
  geom_line(aes(position, WT), color = "black")



his_dt <- readRDS("doc/his_dt.rds")
residueMean <- his_dt[WT1_A_rpc >= 1 & WT2_A_rpc >= 1 & 
                       AT1_A_rpc >= 1 & AT2_A_rpc >= 1 & 
                       WT1_A_sum >= 64 & WT2_A_sum >= 64 & 
                       AT1_A_sum >= 64 & AT2_A_sum >= 64 & 
                       position > 20 & stopdist < -20]
aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')

WT1 <- NULL
WT2 <- NULL
AT1 <- NULL
AT2 <- NULL
for (i in aa) {
  WT1a <- mean(residueMean[residue == i]$WT1_A_pause)
  WT1 <- c(WT1, WT1a)
  WT2a <- mean(residueMean[residue == i]$WT2_A_pause)
  WT2 <- c(WT2, WT2a)
  AT1a <- mean(residueMean[residue == i]$AT1_A_pause)
  AT1 <- c(AT1, AT1a)
  AT2a <- mean(residueMean[residue == i]$AT2_A_pause)
  AT2 <- c(AT2, AT2a)
}
pauseMean <- data.table(aa, WT1, WT2, AT1, AT2)
pauseMean[, WT := (WT1 + WT2) / 2]
pauseMean[, AT := (AT1 + AT2) / 2]

plot <- ggplot(pauseMean, aes(WT1, WT2, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) + 
  xlim(0.4,1.75) + ylim(0.4,1.75)
plot <- plot + theme_classic(20) + labs(y = "WT2 pause score", x = "WT1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Green.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMean, aes(AT1, AT2, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5)
plot <- plot + theme_classic(20) + labs(y = "AT2 pause score", x = "AT1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Green_AT.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMean, aes(WT, AT, color = aa)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5)
plot <- plot + theme_classic(20) + labs(y = "AT pause score", x = "WT pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Green_ATvsWT.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
