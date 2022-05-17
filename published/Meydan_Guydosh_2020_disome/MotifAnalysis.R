library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
disome_dt <- readRDS("disome_dt.rds")

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

KR3of3_dt <- disome_dt[KR3of3, allow.cartesian = TRUE]
KR4of4_dt <- disome_dt[KR4of4, allow.cartesian = TRUE]
KR5of5_dt <- disome_dt[KR5of5, allow.cartesian = TRUE]
KR6of6_dt <- disome_dt[KR6of6, allow.cartesian = TRUE]
KR3of3_dt[, adjusted := position - polybasic]
KR4of4_dt[, adjusted := position - polybasic]
KR5of5_dt[, adjusted := position - polybasic]
KR6of6_dt[, adjusted := position - polybasic]


### Normalize across polybasic interval ###
test <- KR6of6_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("name"))
setkeyv(KR6of6_dt, c("name"))
test[, WT_di1_A_rpc_adjusted := mean(WT_di1_A), by = name]
test[, WT_di2_A_rpc_adjusted := mean(WT_di2_A), by = name]
test[, WT_di3_A_rpc_adjusted := mean(WT_di3_A), by = name]
test[, WT_mono1_A_rpc_adjusted := mean(WT_mono1_A), by = name]
test[, WT_mono2_A_rpc_adjusted := mean(WT_mono2_A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(KR6of6_dt$name, test$name))
KR6of6_dt <- cbind(KR6of6_dt, WT_di1_A_rpc_adjusted = test[i]$WT_di1_A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT_di2_A_rpc_adjusted = test[i]$WT_di2_A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT_di3_A_rpc_adjusted = test[i]$WT_di3_A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT_mono1_A_rpc_adjusted = test[i]$WT_mono1_A_rpc_adjusted)
KR6of6_dt <- cbind(KR6of6_dt, WT_mono2_A_rpc_adjusted = test[i]$WT_mono2_A_rpc_adjusted)
KR6of6_dt[, WT_di1_A_norm := WT_di1_A / WT_di1_A_rpc_adjusted]
KR6of6_dt[, WT_di2_A_norm := WT_di2_A / WT_di2_A_rpc_adjusted]
KR6of6_dt[, WT_di3_A_norm := WT_di3_A / WT_di3_A_rpc_adjusted]
KR6of6_dt[, WT_mono1_A_norm := WT_mono1_A / WT_mono1_A_rpc_adjusted]
KR6of6_dt[, WT_mono2_A_norm := WT_mono2_A / WT_mono2_A_rpc_adjusted]
# saveRDS(KR6of6_dt, "KR6of6_dt.rds")


### Plot ribosome occupancy around region of interest ###
KR3of3_dt <- readRDS("KR3of3_dt.rds")
KR4of4_dt <- readRDS("KR4of4_dt.rds")
KR5of5_dt <- readRDS("KR5of5_dt.rds")
KR6of6_dt <- readRDS("KR6of6_dt.rds")

plot <- ggplot(data = KR6of6_dt[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1 & WT_di3_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT6 = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=3, center=T))]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT6), fun.y = "mean", geom = "line", size=1.25, color = 'red') + 
  stat_summary(data = KR5of5_dt[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT5 = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=3, center=T))], 
               aes(adjusted, WT5), fun.y = "mean", geom = "line", size=1.25, color = 'orange') + 
  stat_summary(data = KR4of4_dt[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT4 = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=3, center=T))], 
               aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'blue') + 
  stat_summary(data = KR3of3_dt[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT3 = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=3, center=T))], 
               aes(adjusted, WT3), fun.y = "mean", geom = "line", size=1.25, color = 'black')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/disome_polybasic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = KR6of6_dt[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT6 = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=3, center=T))]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT6), fun.y = "mean", geom = "line", size=1.25, color = 'red') + 
  stat_summary(data = KR5of5_dt[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT5 = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=3, center=T))], 
               aes(adjusted, WT5), fun.y = "mean", geom = "line", size=1.25, color = 'orange') + 
  stat_summary(data = KR4of4_dt[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT4 = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=3, center=T))], 
               aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'blue') + 
  stat_summary(data = KR3of3_dt[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                .(adjusted, WT3 = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=3, center=T))], 
               aes(adjusted, WT3), fun.y = "mean", geom = "line", size=1.25, color = 'black')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/monosome_polybasic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


KR3of3_dt_up <- KR3of3_dt[, c(1:48)]
test <- KR3of3_dt_up[adjusted >= -50 & adjusted <= -25]
setkeyv(test, c("name"))
setkeyv(KR3of3_dt_up, c("name"))
test[, WT_di1_A_rpc_adjusted := mean(WT_di1_A), by = name]
test[, WT_di2_A_rpc_adjusted := mean(WT_di2_A), by = name]
test[, WT_di3_A_rpc_adjusted := mean(WT_di3_A), by = name]
test[, WT_mono1_A_rpc_adjusted := mean(WT_mono1_A), by = name]
test[, WT_mono2_A_rpc_adjusted := mean(WT_mono2_A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(KR3of3_dt_up$name, test$name))
KR3of3_dt_up <- cbind(KR3of3_dt_up, WT_di1_A_rpc_adjusted = test[i]$WT_di1_A_rpc_adjusted)
KR3of3_dt_up <- cbind(KR3of3_dt_up, WT_di2_A_rpc_adjusted = test[i]$WT_di2_A_rpc_adjusted)
KR3of3_dt_up <- cbind(KR3of3_dt_up, WT_di3_A_rpc_adjusted = test[i]$WT_di3_A_rpc_adjusted)
KR3of3_dt_up <- cbind(KR3of3_dt_up, WT_mono1_A_rpc_adjusted = test[i]$WT_mono1_A_rpc_adjusted)
KR3of3_dt_up <- cbind(KR3of3_dt_up, WT_mono2_A_rpc_adjusted = test[i]$WT_mono2_A_rpc_adjusted)
KR3of3_dt_up[, WT_di1_A_norm := WT_di1_A / WT_di1_A_rpc_adjusted]
KR3of3_dt_up[, WT_di2_A_norm := WT_di2_A / WT_di2_A_rpc_adjusted]
KR3of3_dt_up[, WT_di3_A_norm := WT_di3_A / WT_di3_A_rpc_adjusted]
KR3of3_dt_up[, WT_mono1_A_norm := WT_mono1_A / WT_mono1_A_rpc_adjusted]
KR3of3_dt_up[, WT_mono2_A_norm := WT_mono2_A / WT_mono2_A_rpc_adjusted]

plot <- ggplot(data = KR6of6_dt_up[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1 & WT_di3_A_rpc_adjusted >= 0.1,
                                   .(adjusted, WT6 = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=3, center=T))]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT6), fun.y = "mean", geom = "line", size=1.25, color = 'red') + 
  stat_summary(data = KR5of5_dt_up[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1,
                                   .(adjusted, WT5 = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=3, center=T))], 
               aes(adjusted, WT5), fun.y = "mean", geom = "line", size=1.25, color = 'orange') +
  stat_summary(data = KR4of4_dt_up[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1,
                                   .(adjusted, WT4 = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=3, center=T))], 
               aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'blue') + 
  stat_summary(data = KR3of3_dt_up[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1,
                                   .(adjusted, WT3 = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=3, center=T))], 
               aes(adjusted, WT3), fun.y = "mean", geom = "line", size=1.25, color = 'black')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/disome_polybasic_upstreamNorm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = KR6of6_dt_up[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                   .(adjusted, WT6 = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=3, center=T))]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT6), fun.y = "mean", geom = "line", size=1.25, color = 'red') + 
  stat_summary(data = KR5of5_dt_up[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                   .(adjusted, WT5 = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=3, center=T))], 
               aes(adjusted, WT5), fun.y = "mean", geom = "line", size=1.25, color = 'orange') + 
  stat_summary(data = KR4of4_dt_up[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                   .(adjusted, WT4 = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=3, center=T))], 
               aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'blue') + 
  stat_summary(data = KR3of3_dt_up[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                   .(adjusted, WT3 = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=3, center=T))], 
               aes(adjusted, WT3), fun.y = "mean", geom = "line", size=1.25, color = 'black')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/monosome_polybasic_upstreamNorm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### tripeptide motifs ###
library(tictoc)

# temp <- disome_dt[WT_mono1_A_rpc >= 0.5 & WT_mono2_A_rpc >= 0.5 &
#                     WT_di1_A_rpc >= 0.5 & WT_di2_A_rpc >= 0.5 & WT_di3_A_rpc >= 0.5 &
#                     WT_mono1_A_sum >= 64 & WT_mono2_A_sum >= 64 &
#                     WT_di1_A_sum >= 64 & WT_di2_A_sum >= 64 & WT_di3_A_sum >= 64 &
#                  WT_mono1_A > 0 & !is.na(WT_mono1_A_pause) & WT_mono2_A > 0 & !is.na(WT_mono2_A_pause) &
#                  WT_di1_A > 0 & !is.na(WT_di1_A_pause) & WT_di2_A > 0 & !is.na(WT_di2_A_pause) & WT_di3_A > 0 & !is.na(WT_di3_A_pause) &
#                  position > 20 & stopdist < -20]
temp <- disome_dt[WT_mono1_A > 0 & !is.na(WT_mono1_A_pause) & WT_mono2_A > 0 & !is.na(WT_mono2_A_pause) &
                    WT_di1_A > 0 & !is.na(WT_di1_A_pause) & WT_di2_A > 0 & !is.na(WT_di2_A_pause) & WT_di3_A > 0 & !is.na(WT_di3_A_pause) &
                    position > 20 & stopdist < -20]
motifs <- unique(temp[!is.na(motif3)]$motif3)

tic("runtime")
motif <- NULL
count <- NULL
mono <- NULL
di <- NULL
p_di_mono <- NULL
for (i in motifs) {
  print(match(i, motifs))
  motif <- c(motif, i)
  count1 <- length(which(temp$motif3 == i))
  count <- c(count, count1)
  mono.1 <- mean(temp[motif3 == i]$WT_mono)
  mono <- c(mono, mono.1)
  di.1 <- mean(temp[motif3 == i]$WT_di)
  di <- c(di, di.1)
  test_di <- wilcox.test(temp[motif3 == i]$WT_mono, temp[motif3 == i]$WT_di)
  p_di_mono.0 <- test_di$p.value
  p_di_mono <- c(p_di_mono, p_di_mono.0)
}
toc()
motif_dt <- data.table(motif = motif, count = count, 
                       mono_pause = mono, di_pause = di,
                       p_di_mono = p_di_mono)
motif_dt[, ratio := di_pause / mono_pause]
motif_dt[, padj := p.adjust(p_di_mono, method = "BH")]

ggplot(motif_dt[count >= 5], aes(mono_pause, di_pause)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) + coord_cartesian(ylim = c(0,30), xlim = c(0,30))
View(motif_dt)

library(Logolas)
bg <- summary(disome_dt[position > 20 & stopdist < -20]$residue) / 
  length(disome_dt[position > 20 & stopdist < -20]$residue)
bg <- bg[c(1:19,21)]
logomaker(motif_dt[count >= 5 & ratio > 1]$motif, type = "EDLogo", bg = bg, color_seed = 6)


temp <- disome_dt[motif3 == "PPP", c(1:2,6)]
PPP_dt <- disome_dt[temp, allow.cartesian=TRUE]
PPP_dt[, adjusted := position - i.position]
test <- PPP_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("i.ID"))
setkeyv(PPP_dt, c("i.ID"))
test[, WT_di1_A_rpc_adjusted := mean(WT_di1_A), by = i.ID]
test[, WT_di2_A_rpc_adjusted := mean(WT_di2_A), by = i.ID]
test[, WT_di3_A_rpc_adjusted := mean(WT_di3_A), by = i.ID]
test[, WT_mono1_A_rpc_adjusted := mean(WT_mono1_A), by = i.ID]
test[, WT_mono2_A_rpc_adjusted := mean(WT_mono2_A), by = i.ID]
test <- test[, .SD[which.min(position)], by = i.ID]
i <- cbind(match(PPP_dt$i.ID, test$i.ID))
PPP_dt <- cbind(PPP_dt, WT_di1_A_rpc_adjusted = test[i]$WT_di1_A_rpc_adjusted)
PPP_dt <- cbind(PPP_dt, WT_di2_A_rpc_adjusted = test[i]$WT_di2_A_rpc_adjusted)
PPP_dt <- cbind(PPP_dt, WT_di3_A_rpc_adjusted = test[i]$WT_di3_A_rpc_adjusted)
PPP_dt <- cbind(PPP_dt, WT_mono1_A_rpc_adjusted = test[i]$WT_mono1_A_rpc_adjusted)
PPP_dt <- cbind(PPP_dt, WT_mono2_A_rpc_adjusted = test[i]$WT_mono2_A_rpc_adjusted)
PPP_dt[, WT_di1_A_norm := WT_di1_A / WT_di1_A_rpc_adjusted]
PPP_dt[, WT_di2_A_norm := WT_di2_A / WT_di2_A_rpc_adjusted]
PPP_dt[, WT_di3_A_norm := WT_di3_A / WT_di3_A_rpc_adjusted]
PPP_dt[, WT_mono1_A_norm := WT_mono1_A / WT_mono1_A_rpc_adjusted]
PPP_dt[, WT_mono2_A_norm := WT_mono2_A / WT_mono2_A_rpc_adjusted]

plot <- ggplot(data = PPP_dt[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1 & WT_di3_A_rpc_adjusted >= 0.1,
                                .(adjusted, di = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=1, center=T))]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, di), fun.y = "mean", geom = "line", size=1.25, color = 'red') + 
  stat_summary(data = PPP_dt[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                                .(adjusted, mono = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=1, center=T))],
               aes(adjusted, mono), fun.y = "mean", geom = "line", size=1.25, color = 'gray30')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/disome_PPP.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


disome_dt[, motif2 := substring(motif3,2,3)]
temp <- disome_dt[motif2 == "PP", c(1:2,6)]
PP_dt <- disome_dt[temp, allow.cartesian=TRUE]
PP_dt[, adjusted := position - i.position]
test <- PP_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("i.ID"))
setkeyv(PP_dt, c("i.ID"))
test[, WT_di1_A_rpc_adjusted := mean(WT_di1_A), by = i.ID]
test[, WT_di2_A_rpc_adjusted := mean(WT_di2_A), by = i.ID]
test[, WT_di3_A_rpc_adjusted := mean(WT_di3_A), by = i.ID]
test[, WT_mono1_A_rpc_adjusted := mean(WT_mono1_A), by = i.ID]
test[, WT_mono2_A_rpc_adjusted := mean(WT_mono2_A), by = i.ID]
test <- test[, .SD[which.min(position)], by = i.ID]
i <- cbind(match(PP_dt$i.ID, test$i.ID))
PP_dt <- cbind(PP_dt, WT_di1_A_rpc_adjusted = test[i]$WT_di1_A_rpc_adjusted)
PP_dt <- cbind(PP_dt, WT_di2_A_rpc_adjusted = test[i]$WT_di2_A_rpc_adjusted)
PP_dt <- cbind(PP_dt, WT_di3_A_rpc_adjusted = test[i]$WT_di3_A_rpc_adjusted)
PP_dt <- cbind(PP_dt, WT_mono1_A_rpc_adjusted = test[i]$WT_mono1_A_rpc_adjusted)
PP_dt <- cbind(PP_dt, WT_mono2_A_rpc_adjusted = test[i]$WT_mono2_A_rpc_adjusted)
PP_dt[, WT_di1_A_norm := WT_di1_A / WT_di1_A_rpc_adjusted]
PP_dt[, WT_di2_A_norm := WT_di2_A / WT_di2_A_rpc_adjusted]
PP_dt[, WT_di3_A_norm := WT_di3_A / WT_di3_A_rpc_adjusted]
PP_dt[, WT_mono1_A_norm := WT_mono1_A / WT_mono1_A_rpc_adjusted]
PP_dt[, WT_mono2_A_norm := WT_mono2_A / WT_mono2_A_rpc_adjusted]

plot <- ggplot(data = PP_dt[WT_di1_A_rpc_adjusted >= 0.1 & WT_di2_A_rpc_adjusted >= 0.1 & WT_di3_A_rpc_adjusted >= 0.1,
                    .(adjusted, di = movingAverage(((WT_di1_A_norm + WT_di2_A_norm + WT_di3_A_norm) / 3), n=1, center=T))]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, di), fun.y = "mean", geom = "line", size=1.25, color = 'red') + 
  stat_summary(data = PP_dt[WT_mono1_A_rpc_adjusted >= 0.1 & WT_mono2_A_rpc_adjusted >= 0.1,
                            .(adjusted, mono = movingAverage(((WT_mono1_A_norm + WT_mono2_A_norm) / 2), n=1, center=T))],
               aes(adjusted, mono), fun.y = "mean", geom = "line", size=1.25, color = 'gray30')
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/disome_PP.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

