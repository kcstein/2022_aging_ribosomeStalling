library(data.table)
library(ggplot2)
source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
GH_dt <- readRDS("GH_dt.rds")
GHA_dt <- readRDS("GHA_dt.rds")


ggplot(aes(position, occupancy), data=GH_dt[orf == "YJR009C", .(position, occupancy = movingAverage(GH0_pause, n=1, center=T),
                                                                    occupancy2 = movingAverage(GH4_pause, n=1, center=T))]) + 
  geom_vline(xintercept = 310, color = "gray50", linetype = "dashed") +
  geom_line(aes(position, occupancy2), color = 'red', size = 1.25) +
  geom_line(color = "gray30", size = 1.25)
  
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (pause score)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/GHreporter_pause.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(aes(position, occupancy), data=GH_dt[orf == "GFP_HIS3", .(position, occupancy = movingAverage(((GH0_1A_rpm + GH0_2A_rpm) / 2), n=5, center=T),
                                                                       occupancy2 = movingAverage(((GH4_1A_rpm + GH4_2A_rpm) / 2), n=5, center=T))]) + 
  #geom_vline(xintercept = 241, color = "gray50", linetype = "dashed") +
  geom_line(color = "gray40", size = 1.25) +
  geom_line(aes(position, occupancy2), color = '#E7A427', size = 1.25)
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (rpm)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/GHreporter_rpm.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(aes(position, occupancy), data=GHA_dt[orf == "GFP_AAA12_HIS3", .(position, occupancy = movingAverage(GHA0_pause, n=5, center=T),
                                                                       occupancy2 = movingAverage(GHA4_pause, n=5, center=T))]) + 
  #geom_vline(xintercept = 241, color = "gray50", linetype = "dashed") +
  geom_line(color = "gray40", size = 1.25) +
  geom_line(aes(position, occupancy2), color = '#E7A427', size = 1.25)
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (pause score)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/K12reporter_pause.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=GHA_dt[orf == "GFP_AAA12_HIS3", .(position, occupancy = movingAverage(((GHA0_1A_rpm + GHA0_2A_rpm) / 2), n=5, center=T),
                                                                       occupancy2 = movingAverage(((GHA4_1A_rpm + GHA4_2A_rpm) / 2), n=5, center=T))]) + 
  #geom_vline(xintercept = 241, color = "gray50", linetype = "dashed") +
  geom_line(color = "gray40", size = 1.25) +
  geom_line(aes(position, occupancy2), color = '#E7A427', size = 1.25)
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/K12reporter_rpm.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(aes(position, occupancy), data=GHA_dt[orf == "GFP_AAA12_HIS3", .(position, occupancy = movingAverage(GHA0_pause, n=1, center=T),
                                                                       occupancy2 = movingAverage(GHA4_pause, n=1, center=T))]) + 
  geom_vline(xintercept = 241, color = "gray50", linetype = "dashed") +
  geom_line(color = "gray30", size = 1.25) +
  geom_line(aes(position, occupancy2), color = 'red', size = 1.25) +
  coord_cartesian(xlim = c(221,261))
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (pause score)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/K12reporter_pause_insert.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=GHA_dt[orf == "GFP_AAA12_HIS3", .(position, occupancy = movingAverage(((GHA0_1A_rpm + GHA0_2A_rpm) / 2), n=1, center=T),
                                                                               occupancy2 = movingAverage(((GHA4_1A_rpm + GHA4_2A_rpm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 241, color = "gray50", linetype = "dashed") +
  geom_line(color = "gray30", size = 1.25) +
  geom_line(aes(position, occupancy2), color = 'red', size = 1.25) +
  coord_cartesian(xlim = c(221,261))
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/K12reporter_rpm_insert.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

