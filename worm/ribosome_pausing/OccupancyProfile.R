library(ggplot2)


### Start and stop codon analysis ###
# Plot traces from start codon
N2_dtA <- readRDS("N2_dtA.rds")
ggplot(data = N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & 
                       D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & length >= 300,
                     .(position, D1 = ((D1_1A_pause + D1_2A_pause) / 2),
                       D12 = ((D12_1A_pause + D12_2A_pause) / 2))]) + xlim(-7, 300) +
  stat_summary(aes(position, D1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(position, D12), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(position, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(position, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'red')

# Plot from the stop codon
ggplot(data = N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & 
                       D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & length >= 300,
                     .(stopdist, D1 = ((D1_1A_pause + D1_2A_pause) / 2),
                       D12 = ((D12_1A_pause + D12_2A_pause) / 2))]) + xlim(-50, 7) +
  stat_summary(aes(stopdist, D1), fun.y = "mean", geom = "line", size=0.5, color = 'black') +
  stat_summary(aes(stopdist, D12), fun.y = "mean", geom = "line", size=0.5, color = 'red') +
  stat_summary(aes(stopdist, D1), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(aes(stopdist, D12), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'red')


### Occupancy in single gene ###
stalling_peaks_ageDep <- readRDS("stalling_peaks_ageDep.rds")
plot <- ggplot(aes(position, occupancy), data=N2_dtA[gene == "mars-1", .(position, occupancy = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=5, center=T),
                                                                         occupancy2 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=5, center=T))]) + 
  geom_vline(xintercept = stalling_peaks_ageDep[gene == "mars-1"]$position, color = "gray50", linetype = "dashed") +
  geom_line(aes(position, occupancy2), color = 'red', size = 1.25) +
  geom_line(color = "gray30", size = 1.25)
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/cct-8.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=N2_dtA[gene == "wars-1", .(position, D1 = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=5, center=T), 
                                               D12 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=5, center=T))]) + 
  geom_line(aes(position, D12), color = "#E31A1C", size = 1.25) +
  geom_line(aes(position, D1), color = "gray30", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(53,93))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/wars-1_aa166.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

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
plot <- ggplot(aes(position, occupancy), data=bukau_dt["YJR121W", .(position, occupancy = movingAverage(B_ssb1_I_rpm, n=3, center=T))]) + 
  geom_vline(xintercept = bukau_dataset_dt[orf == "YJR121W"]$peak, color = "blue", linetype = "dashed") +
  geom_line(aes(color = "1"),size = 1.25) + 
  geom_line(aes(position, occupancy2, color = "2"), data=bukau_dt["YJR121W", .(position, occupancy2 = movingAverage(B_ssb1_T_rpm, n=3, center=T))], size = 1.25) +
  scale_color_manual(values = c("1" = "red", "2" = "black"), labels = c("Ssb IP", "Total"), name = "")
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("ATP2")
ggsave("/Users/KevinStein/Desktop/ATP2.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


### Save profiles of set of genes ###
temp <- unique(KR5of5_dt[D1_1A_rpc_adjusted >= 0.1 & D1_2A_rpc_adjusted >= 0.1 & 
                           D12_1A_rpc_adjusted >= 0.1 & D12_2A_rpc_adjusted >= 0.1 & polybasic > 25 & adjusted == 0, 1])
temp1 <- KR5of5_dt[adjusted == 0]
for (i in temp$orf) {
  filename <- paste("/Users/KevinStein/Desktop/KR5of5/", i, ".pdf", sep = "")
  G <- ggplot(aes(position, occupancy), data=N2_dtA[orf == as.character(i), .(position, occupancy = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=5, center=T))]) + 
    geom_vline(data = temp1[orf == as.character(i)], xintercept = temp1[orf == as.character(i)]$polybasic) + 
    geom_line(color = 'black') +
    geom_line(aes(position, occupancy2), data=N2_dtA[orf == as.character(i), .(position, occupancy2 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=5, center=T))], color = 'red')
  ggsave(filename, G, width = 6, height = 4, dpi = 300)
}

temp <- c("Q86D21", "P91128", "P34339", "Q09996", "P29691", "Q18678", "Q22918", 
          "P48166", "Q9NAP9", "Q9U1Q4", "O62431", "P51403", "O01974", "Q10039",
          "Q19722", "Q19825", "Q7JMF0", "O01541", "Q03577", "Q20970", "G5EDY2", "Q9XVE9")
for (i in temp) {
  filename <- paste("/Users/KevinStein/Desktop/translation/", i, ".pdf", sep = "")
  G <- ggplot(aes(position, occupancy), data=N2_dtA[uniprot == i, .(position, occupancy = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=1, center=T))]) + 
    geom_vline(data = stalling_peaks[uniprot == i], xintercept = stalling_peaks[uniprot == i]$peak) + 
    geom_line(color = 'black') +
    geom_line(aes(position, occupancy2), data=N2_dtA[uniprot == i, .(position, occupancy2 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=1, center=T))], color = 'red')
  ggsave(filename, G, width = 6, height = 4, dpi = 300)
}

temp <- unique(stalling_peaks_ageDep[P == "W", 1])
temp <- unique(stalling_peaks_ageDep[residue == "P", 1])
temp1 <- stalling_peaks_ageDep_dt[residue == "P" & adjusted == 0]
for (i in temp$orf) {
  filename <- paste("/Users/KevinStein/Desktop/P/", i, ".pdf", sep = "")
  G <- ggplot(aes(position, occupancy), data=N2_dtA[orf == as.character(i), .(position, occupancy = movingAverage(((D1_1A_pause + D1_2A_pause) / 2), n=5, center=T))]) + 
    geom_vline(data = temp1[orf == as.character(i)], xintercept = temp1[orf == as.character(i)]$peak) + 
    geom_line(color = 'black') +
    geom_line(aes(position, occupancy2), data=N2_dtA[orf == as.character(i), .(position, occupancy2 = movingAverage(((D12_1A_pause + D12_2A_pause) / 2), n=5, center=T))], color = 'red')
  ggsave(filename, G, width = 6, height = 4, dpi = 300)
}
