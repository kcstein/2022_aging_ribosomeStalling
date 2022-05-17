library(ggplot2)
sc_dtA <- readRDS("doc/sc_dtA.rds")


### Start and stop codon analysis ###
# Plot traces from start codon
ggplot(data = sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 & 
                       WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 &
                       sch9_D0_1A_rpc >= 0.5 & sch9_D0_2A_rpc >= 0.5 & 
                       sch9_D4_1A_rpc >= 0.5 & sch9_D4_2A_rpc >= 0.5 & length >= 300,
                     .(position, WT0 = ((WT_D0_1A_pause + WT_D0_2A_pause) / 2),
                       WT4 = ((WT_D4_1A_pause + WT_D4_2A_pause) / 2),
                       S0 = ((sch9_D0_1A_pause + sch9_D0_2A_pause) / 2),
                       S4 = ((sch9_D4_1A_pause + sch9_D4_2A_pause) / 2))]) + xlim(-7, 300) +
  stat_summary(aes(position, S0), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(aes(position, S4), fun.y = "mean", geom = "line", size=0.5, color = 'orange') +
  stat_summary(aes(position, WT0), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(position, WT4), fun.y = "mean", geom = "line", size=0.5, color = 'red')
#stat_summary(aes(position, WT0), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
#     fun.args=list(conf.int=0.5), fill = 'black') +
#stat_summary(aes(position, WT4), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
#            fun.args=list(conf.int=0.5), fill = 'red')


# Plot from the stop codon
ggplot(data = sc_dtA[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 & 
                       WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 &
                       sch9_D0_1A_rpc >= 0.5 & sch9_D0_2A_rpc >= 0.5 & 
                       sch9_D4_1A_rpc >= 0.5 & sch9_D4_2A_rpc >= 0.5 & length >= 300,
                     .(stopdist, WT0 = ((WT_D0_1A_pause + WT_D0_2A_pause) / 2),
                       WT4 = ((WT_D4_1A_pause + WT_D4_2A_pause) / 2),
                       S0 = ((sch9_D0_1A_pause + sch9_D0_2A_pause) / 2),
                       S4 = ((sch9_D4_1A_pause + sch9_D4_2A_pause) / 2))]) + xlim(-50, 7) +
  stat_summary(aes(stopdist, S0), fun.y = "mean", geom = "line", size=0.5, color = 'blue') +
  stat_summary(aes(stopdist, S4), fun.y = "mean", geom = "line", size=0.5, color = 'orange') +
  stat_summary(aes(stopdist, WT0), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(stopdist, WT4), fun.y = "mean", geom = "line", size=0.5, color = 'red')
#stat_summary(aes(stopdist, WT0), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
#    fun.args=list(conf.int=0.5), fill = 'black') +
#stat_summary(aes(stopdist, WT4), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
#           fun.args=list(conf.int=0.5), fill = 'red')

ggplot(data = sc_dt[WT_D0_1A_rpc >= 0.5 & WT_D0_2A_rpc >= 0.5 & 
                       WT_D4_1A_rpc >= 0.5 & WT_D4_2A_rpc >= 0.5 & length >= 100,
                     .(stopdist, WT0 = ((WT_D0_1A_pause + WT_D0_2A_pause) / 2),
                       WT4 = ((WT_D4_1A_pause + WT_D4_2A_pause) / 2))]) + xlim(-50, 195) +
  stat_summary(aes(stopdist, WT0), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(stopdist, WT4), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  coord_cartesian(xlim = c(-20,20))

utr3 <- sc_dt[stopdist > 0 & stopdist < 50]
utr3[, utr3 := sum(WT_D0_1A)]

### Occupancy in single gene ###
ggplot(aes(position, occupancy), data=expression1["YGL116W", .(position, occupancy = movingAverage(ribo_rpm, n=5, center=T),
                                                           occupancy2 = movingAverage(tric_rpm, n=5, center=T),
                                                           occupancy3 = movingAverage(ribo_tpm, n=5, center=T),
                                                           occupancy4 = movingAverage(tric_tpm, n=5, center=T))]) +
  geom_line(aes(position, occupancy2), color = 'red') + geom_line() + geom_line(aes(position, occupancy3), color = 'blue') + geom_line(aes(position, occupancy4), color = 'orange') + coord_cartesian(ylim = c(-1,2))

# Padj
plot <- ggplot(aes(position, tric_ribo), data=tric_fishers["YJR121W", .(position, tric_ribo = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_line(color = "gray50") +
  geom_point(data = tric_fishers["YJR121W" & (tric_sum > (1.5*tric_rpc)), .(position, tric_ribo = movingAverage(tric_odds, n=5, center=T))], aes(color = tric_fishers[orf == "YJR121W"]$tric_padj)) + 
  scale_color_distiller(type = "seq", palette = "Reds", limits = c(0,0.05), name = "Adjusted p-value", guide = guide_colorbar(barwidth = 0.75, barheight = 3)) +
  scale_x_continuous(breaks = c(0,100,200,300,400,500)) + 
  scale_y_continuous(limits = c(-10,31), breaks = c(0,10,20,30))
scale_color_gradient(name = "Adjusted p-value", low = "#67000D", high = "#FFF5F0", na.value = "gray50", limits = c(0,1), guide = guide_colorbar(barwidth = 0.75, barheight = 3))

# Two different datasets
ggplot(aes(position, occupancy), data=tric_dt[orf == "YNL317W", .(position, occupancy = movingAverage(tric_rpm, n=5, center=T))]) + geom_line(color = 'red') +
  geom_line(aes(position, occupancy2), data=tric_dt[orf == "YNL317W", .(position, occupancy2 = movingAverage(ribo_rpm, n=5, center=T))], color = 'black')

# Pretty
ggplot(aes(position, occupancy), data=sc_dtA[orf == "YDR190C", .(position, occupancy = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=1, center=T),
                                                                         occupancy2 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 220, color = "gray50", linetype = "dashed") +
  geom_vline(xintercept = 460, color = "gray50", linetype = "dashed") +
  #geom_vline(xintercept = 85, color = "gray50", linetype = "dashed") +
  #geom_vline(xintercept = WT_stalling_peaks_ageDep[orf == "YJL008C"]$position, color = "gray50", linetype = "dashed") +
  geom_line(aes(position, occupancy2), color = 'red', size = 1.25) +
  geom_line(color = "gray30", size = 1.25)
  
ggsave("/Users/KevinStein/Desktop/CCT8.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
  #geom_line(aes(color = "1"),size = 1.25) + 
  #geom_line(aes(position, occupancy2, color = "2"), data=bukau_dt["YJR121W", .(position, occupancy2 = movingAverage(B_ssb1_T_rpm, n=3, center=T))], size = 1.25) +
  #scale_color_manual(values = c("1" = "red", "2" = "black"), labels = c("Ssb IP", "Total"), name = "")
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("RNQ1")
ggsave("/Users/KevinStein/Desktop/RNQ1.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


### Save profiles of set of genes ###
KR6of6_dt <- readRDS("doc/KR6of6_dt.rds")
temp <- unique(KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                           WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                           polybasic > 25 & (polybasic - length) < -25 & adjusted == 0, 1])
temp1 <- KR6of6_dt[adjusted == 0]
for (i in temp$orf) {
  filename <- paste("/Users/KevinStein/Desktop/KR6of6_yeast/", i, ".pdf", sep = "")
  G <- ggplot(aes(position, occupancy), data=sc_dtA[orf == as.character(i), .(position, occupancy = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=5, center=T))]) + 
    geom_vline(data = temp1[orf == as.character(i)], xintercept = temp1[orf == as.character(i)]$polybasic) + 
    geom_line(color = 'black') +
    geom_line(aes(position, occupancy2), data=sc_dtA[orf == as.character(i), .(position, occupancy2 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=5, center=T))], color = 'red')
  ggsave(filename, G, width = 6, height = 4, dpi = 300)
}

WT_stalling_peaks_ageDep <- readRDS("WT_stalling_peaks_ageDep.rds")
WT_stalling_peaks_ageDep_dt <- readRDS("WT_stalling_peaks_ageDep_dt.rds")
temp <- unique(WT_stalling_peaks_ageDep[residue == "Y", 1])
temp1 <- WT_stalling_peaks_ageDep_dt[residue == "Y" & adjusted == 0]
for (i in temp$orf) {
  filename <- paste("/Users/KevinStein/Desktop/Yyeast/", i, ".pdf", sep = "")
  G <- ggplot(aes(position, occupancy), data=sc_dtA[orf == as.character(i), .(position, occupancy = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=5, center=T))]) + 
    geom_vline(data = temp1[orf == as.character(i)], xintercept = temp1[orf == as.character(i)]$peak) + 
    geom_line(color = 'black') +
    geom_line(aes(position, occupancy2), data=sc_dtA[orf == as.character(i), .(position, occupancy2 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=5, center=T))], color = 'red')
  ggsave(filename, G, width = 6, height = 4, dpi = 300)
}

