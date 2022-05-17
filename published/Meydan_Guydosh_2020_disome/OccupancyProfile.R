library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
disome_dt <- readRDS("disome_dt.rds")


ggplot(aes(position, occupancy), data=disome_dt[orf == "YEL056W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=1, center=T),
                                                                 occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=1, center=T))]) + 
  #geom_vline(xintercept = 227, color = "gray50", linetype = "dashed") +
  geom_line(aes(position, occupancy2), color = 'red', size = 1.25) +
  geom_line(color = "gray30", size = 1.25)
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("RNQ1")
ggsave("/Users/KevinStein/Desktop/RNQ1.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


ggplot(aes(position, occupancy), data=disome_dt[orf == "YOR202W", .(position, occupancy = movingAverage(WT_His3_di_A_rpm, n=1, center=T))]) + 
  #geom_vline(xintercept = 227, color = "gray50", linetype = "dashed") +
  geom_line(color = "gray30", size = 1.25)

ggplot(aes(position, occupancy), data=disome_dt[orf == "YOR202W", .(position, occupancy = movingAverage(WT_His3_di_A_rpm, n=1, center=T),
                                                                    occupancy2 = movingAverage(WT_His3_tri_A_rpm, n=1, center=T))]) + 
  #geom_vline(xintercept = 227, color = "gray50", linetype = "dashed") +
  geom_line(aes(position, occupancy2), color = 'red', size = 1.25) +
  geom_line(color = "gray30", size = 1.25)


