### Align ribosome occupancy at start of polybasic stretch ###
dir <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/"
files <- list.files(dir, pattern = "*_only.csv")
samples <- sub("\\_only.csv", '', files)

for(i in 1:length(files)) {
  print(files[i])
  j <- as.data.table(lapply(paste0(dir, files[i]), read.csv, header = TRUE, stringsAsFactors = TRUE))
  setnames(j, c("orf", "polybasic", "name"))
  setkeyv(j, c("orf"))
  assign(samples[i], j)
}

KR3of3_dt <- N2_dtA[KR3of3, allow.cartesian=TRUE]
KR4of4_dt <- N2_dtA[KR4of4, allow.cartesian=TRUE]
KR5of5_dt <- N2_dtA[KR5of5, allow.cartesian=TRUE]
KR6of6_dt <- N2_dtA[KR6of6, allow.cartesian=TRUE]
KR3of3_dt[, adjusted := position - polybasic]
KR4of4_dt[, adjusted := position - polybasic]
KR5of5_dt[, adjusted := position - polybasic]
KR6of6_dt[, adjusted := position - polybasic]

### Normalize across polybasic interval ###
test <- KR6of6_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("name"))
setkeyv(KR6of6_dt, c("name"))
test[, D1_1A_rpc_adjusted := mean(D1_1A), by = name]
test[, D1_2A_rpc_adjusted := mean(D1_2A), by = name]
test[, D12_1A_rpc_adjusted := mean(D12_1A), by = name]
test[, D12_2A_rpc_adjusted := mean(D12_2A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(KR6of6_dt$name, test$name))
KR6of6_dt <- cbind(KR6of6_dt, D1_1A_rpc_adjusted = test[i]$D1_1A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, D1_2A_rpc_adjusted = test[i]$D1_2A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, D12_1A_rpc_adjusted = test[i]$D12_1A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, D12_2A_rpc_adjusted = test[i]$D12_2A_rpc_adjusted)
KR6of6_dt[, D1_1A_norm := D1_1A / D1_1A_rpc_adjusted]
KR6of6_dt[, D1_2A_norm := D1_2A / D1_2A_rpc_adjusted]
KR6of6_dt[, D12_1A_norm := D12_1A / D12_1A_rpc_adjusted]
KR6of6_dt[, D12_2A_norm := D12_2A / D12_2A_rpc_adjusted]
# saveRDS(KR6of6_dt, "KR6of6_dt.rds")


### Plot ribosome occupancy around site of interest ###
KR3of3_dt <- readRDS("doc/KR3of3_dt.rds")
KR4of4_dt <- readRDS("doc/KR4of4_dt.rds")
KR5of5_dt <- readRDS("doc/KR5of5_dt.rds")
KR6of6_dt <- readRDS("doc/KR6of6_dt.rds")

# average replicates, adjusted rpc0.1
ggplot(data = KR4of4_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & polybasic > 25,
                        .(adjusted, D1A_pause = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=5, center=T),
                          D12A_pause = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=5, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1A_pause), fun.y = "median", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D12A_pause), fun.y = "median", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(adjusted, D1A_pause), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(adjusted, D12A_pause), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'red')

# final
ggplot(data = KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                          D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1,
                        .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=5, center=T),
                          D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=5, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'red')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/Worm_KR3of3_original.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                        D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 70 & (polybasic - length) < -25 & adjusted == 3]$D1_1A_norm,
            KR6of6_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                        D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & 
                        polybasic > 70 & (polybasic - length) < -25 & adjusted == 3]$D12_2A_norm, alternative = 't')

ggplot(data = KR3of3_dt[N2_D12_1A_rpc_adjusted >= 0.1 & N2_D12_2A_rpc_adjusted >= 0.1,
                            .(adjusted, D1_2 = movingAverage(N2_D1_2A_norm, n=3, center=T),
                              D12_2 = movingAverage(N2_D12_2A_norm, n=3, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1_2), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(adjusted, D12_2), fun.y = "mean", geom = "line", size=0.5, color = 'red') + coord_cartesian(ylim = c(0.8,1.4))


ggplot(data = daf_KR4of4_dt[daf2_D12_1A_rpc_adjusted >= 0.1 & daf2_D12_2A_rpc_adjusted >= 0.1,
                        .(adjusted, D1_2 = movingAverage(daf2_D1_2A_norm, n=3, center=T),
                          D12_2 = movingAverage(daf2_D12_2A_norm, n=3, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1_2), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(adjusted, D12_2), fun.y = "mean", geom = "line", size=0.5, color = 'red') + coord_cartesian(ylim = c(0.8,1.4))

ggplot(data = daf_KR4of4_dt[daf16_D12_1A_rpc_adjusted >= 0.1 & daf16_D12_2A_rpc_adjusted >= 0.1,
                            .(adjusted, D1_2 = movingAverage(daf16_D1_2A_norm, n=3, center=T),
                              D12_2 = movingAverage(daf16_D12_2A_norm, n=3, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1_2), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(adjusted, D12_2), fun.y = "mean", geom = "line", size=0.5, color = 'red') + coord_cartesian(ylim = c(0.8,1.4))

KR5of5_dt <- readRDS("doc/KR5of5_dt.rds")
ggplot(data = KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                          D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1,
                        .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=1, center=T),
                          D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=1, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=0.5, color = 'red')

ggplot(data = KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                          D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1,
                        N2_D6_1A_rpc_adjusted >= 0.1 & N2_D6_2A_rpc_adjusted >= 0.1,
                        .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=3, center=T),
                          D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=3, center=T),
                          D6_1 = movingAverage(N2_D6_1A_norm, n=3, center=T),
                          D6_2 = movingAverage(N2_D6_2A_norm, n=3, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(adjusted, D6_1), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(aes(adjusted, D6_2), fun.y = "mean", geom = "line", size=0.5, color = 'orange')


ggplot(data = KR4of4_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & N2_D6_2A_rpc_adjusted >= 0.1 &
                          D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & position > 25,
                        .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=3, center=T),
                          D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=3, center=T),
                          D6 = movingAverage(N2_D6_2A_norm, n=3, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +  
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(adjusted, D6), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(data = daf_KR4of4_dt[daf16_D1_2A_rpc_adjusted >= 0.1 & daf16_D6_2A_rpc_adjusted >= 0.1 &
                                      daf16_D12_2A_rpc_adjusted >= 0.1 & position > 25,
                                    .(adjusted, daf16_D1 = movingAverage(daf16_D1_2A_norm, n=3, center=T))],
               aes(adjusted, daf16_D1), fun.y = "mean", geom = "line", size=0.5, color = 'gray50') +
  stat_summary(data = daf_KR4of4_dt[daf16_D1_2A_rpc_adjusted >= 0.1 & daf16_D6_2A_rpc_adjusted >= 0.1 &
                                      daf16_D12_2A_rpc_adjusted >= 0.1 & position > 25,
                                    .(adjusted, daf16_D12 = movingAverage(daf16_D12_2A_norm, n=3, center=T))],
               aes(adjusted, daf16_D12), fun.y = "mean", geom = "line", size=0.5, color = 'orange') +
  stat_summary(data = daf_KR4of4_dt[daf16_D1_2A_rpc_adjusted >= 0.1 & daf16_D6_2A_rpc_adjusted >= 0.1 &
                                      daf16_D12_2A_rpc_adjusted >= 0.1 & position > 25,
                                    .(adjusted, daf16_D6 = movingAverage(daf16_D6_2A_norm, n=3, center=T))],
               aes(adjusted, daf16_D6), fun.y = "mean", geom = "line", size=0.5, color = 'green')


R4of4_dt <- readRDS("doc/Archive_Misc/R4of4_dt.rds")
ggplot(data = R4of4_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                          D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1,
                        .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=4, center=T),
                          D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E31A1C') +
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size = 1.25, color = '#E31A1C') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size = 1.25, color = 'black') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())


### RR
RR <- N2_dtA[motif2 == "RR" & position > 15 & stopdist < -15 &
               D1_1A_rpc >= 1 & D1_2A_rpc >= 1 &
               D12_1A_rpc >= 1 & D12_2A_rpc >= 1]
ggplot(RR, aes("1", D1_pause)) + geom_boxplot() +
  geom_boxplot(aes("12", D12_pause)) + scale_y_log10()
wilcox.test(RR$D1_pause, RR$D12_pause)
RR1 <- RR[, c(1:2,6)]
setnames(RR1, c("orf", "peak", "name"))
RR_dt <- N2_dtA[RR1, allow.cartesian=TRUE]
RR_dt[, adjusted := position - peak]

### Normalize across polybasic interval ###
test <- RR_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("name"))
setkeyv(RR_dt, c("name"))
test[, D1_1A_rpc_adjusted := mean(D1_1A), by = name]
test[, D1_2A_rpc_adjusted := mean(D1_2A), by = name]
test[, D12_1A_rpc_adjusted := mean(D12_1A), by = name]
test[, D12_2A_rpc_adjusted := mean(D12_2A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(RR_dt$name, test$name))
RR_dt <- cbind(RR_dt, D1_1A_rpc_adjusted = test[i]$D1_1A_rpc_adjusted)
RR_dt <- cbind(RR_dt, D1_2A_rpc_adjusted = test[i]$D1_2A_rpc_adjusted)
RR_dt <- cbind(RR_dt, D12_1A_rpc_adjusted = test[i]$D12_1A_rpc_adjusted)
RR_dt <- cbind(RR_dt, D12_2A_rpc_adjusted = test[i]$D12_2A_rpc_adjusted)
RR_dt[, D1_1A_norm := D1_1A / D1_1A_rpc_adjusted]
RR_dt[, D1_2A_norm := D1_2A / D1_2A_rpc_adjusted]
RR_dt[, D12_1A_norm := D12_1A / D12_1A_rpc_adjusted]
RR_dt[, D12_2A_norm := D12_2A / D12_2A_rpc_adjusted]
# saveRDS(RR_dt, "RR_dt.rds")

plot <- ggplot(data = RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                      D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1,
                    .(adjusted, D1 = movingAverage(((D1_1A_norm + D1_2A_norm) / 2), n=3, center=T),
                      D12 = movingAverage(((D12_1A_norm + D12_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(adjusted, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E31A1C') +
  stat_summary(aes(adjusted, D12), fun.y = "mean", geom = "line", size = 1.25, color = '#E31A1C') +
  stat_summary(aes(adjusted, D1), fun.y = "mean", geom = "line", size = 1.25, color = 'black') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
wilcox.test(RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                        D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 & 
                        adjusted == 0]$D1_2A_norm,
            RR_dt[D1_1A_rpc_adjusted >= 1 & D1_2A_rpc_adjusted >= 1 & 
                        D12_1A_rpc_adjusted >= 1 & D12_2A_rpc_adjusted >= 1 & 
                        adjusted == 0]$D12_2A_norm, alternative = 't')

