
### Sis1 KR6of10 ###
plot <- ggplot(data=sc_dtA[orf == "YNL007C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_rect(xmin = 177, xmax = 217, ymin = -Inf, ymax = Inf, fill = "#F2CF8C", alpha = 0.1) +
  geom_vline(xintercept = 197, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expansion())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YNL007C_SIS1_KR6of10aa197.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YNL007C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                       WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 197, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expansion(), limits = c(177,217))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YNL007C_SIS1_KR6of10aa197_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# All polybasic stretches from aged cells in one plot
KR3of3_dt <- readRDS("doc/KR3of3_dt.rds")
KR4of4_dt <- readRDS("doc/KR4of4_dt.rds")
KR5of5_dt <- readRDS("doc/KR5of5_dt.rds")
KR6of6_dt <- readRDS("doc/KR6of6_dt.rds")
plot <- ggplot(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'red') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'red') +
  stat_summary(data = KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'orange') +
  stat_summary(data = KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'orange') +
  stat_summary(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'blue') +
  stat_summary(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'blue') +
  stat_summary(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = 'black') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale()) + coord_cartesian(ylim = c(0.5,2.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/WT_allPolybasic_Day4.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = KR6of6_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'red') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'red') +
  stat_summary(data = KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'orange') +
  stat_summary(data = KR5of5_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'orange') +
  stat_summary(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'blue') +
  stat_summary(data = KR4of4_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'blue') +
  stat_summary(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'black') +
  stat_summary(data = KR3of3_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  polybasic > 25 & (polybasic - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=3, center=T))],
               aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'black') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale()) + coord_cartesian(ylim = c(0.5,2.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from polybasic start (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/WT_allPolybasic_Day0.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



### GCN4 and ribosomal protein ###
sc_dt_UTR <- readRDS("sc_dt_UTR.rds")
setnames(sc_dt_UTR, c("orf","position","codon","length","stopdist","ID","residue",
                         "WT_D0_1A","WT_D0_2A","WT_D2_1A","WT_D2_2A",       "WT_D4_1A",      
"WT_D4_2A",       "WT_D0_1A_sum","WT_D0_1A_rpc","WT_D0_2A_sum",  
"WT_D0_2A_rpc",   "WT_D2_1A_sum","WT_D2_1A_rpc","WT_D2_2A_sum",  
"WT_D2_2A_rpc",   "WT_D4_1A_sum",   "WT_D4_1A_rpc","WT_D4_2A_sum",  
"WT_D4_2A_rpc",   "WT_D0_1A_rpm",   "WT_D0_1A_pause","WT_D0_2A_rpm",  
"WT_D0_2A_pause", "WT_D2_1A_rpm",   "WT_D2_1A_pause","WT_D2_2A_rpm",  
"WT_D2_2A_pause", "WT_D4_1A_rpm",   "WT_D4_1A_pause","WT_D4_2A_rpm",  
"WT_D4_2A_pause"))
plot <- ggplot(data=sc_dtA[orf == "YEL009C", .(position, WT0 = movingAverage(((WT_D0_1A_rpm + WT_D0_2A_rpm) / 2), n=1, center=T), 
                                          WT2 = movingAverage(((WT_D2_1A_rpm + WT_D2_2A_rpm) / 2), n=1, center=T),     
                                          WT4 = movingAverage(((WT_D4_1A_rpm + WT_D4_2A_rpm) / 2), n=1, center=T))]) + 
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT2), color = "#F2CF8C", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expansion())
plot <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (rpm)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YEL009C_GCN4cds.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(data=sc_dt_UTR[orf == "YEL009C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=1, center=T), 
                                          WT2 = movingAverage(((WT_D2_1A_pause + WT_D2_2A_pause) / 2), n=1, center=T),     
                                          WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=1, center=T))]) + 
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  geom_line(aes(position, WT2), color = "#F2CF8C", size = 1.25) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  scale_x_continuous(expand = expansion(), limits = c(-130,0))
  plot <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (rpm)", x = "Codon position") +
    theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
  ggsave("/Users/KevinStein/Desktop/YEL009C_GCN4utr.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
  
  plot <- ggplot(data=sc_dtA[orf == "YOR063W", .(position, WT0 = movingAverage(((WT_D0_1A_rpm + WT_D0_2A_rpm) / 2), n=1, center=T), 
                                                 WT2 = movingAverage(((WT_D2_1A_rpm + WT_D2_2A_rpm) / 2), n=1, center=T),     
                                                 WT4 = movingAverage(((WT_D4_1A_rpm + WT_D4_2A_rpm) / 2), n=1, center=T))]) + 
    geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
    geom_line(aes(position, WT2), color = "#F2CF8C", size = 1.25) +
    geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
    scale_x_continuous(expand = expansion())
  plot <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (rpm)", x = "Codon position") +
    theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
  ggsave("/Users/KevinStein/Desktop/YOR063W_RPL3.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
  

  ### Metagene from start
  sc_dt_UTR[, WT_D0_pause := (WT_D0_1A_pause + WT_D0_2A_pause) / 2]
  sc_dt_UTR[, WT_D2_pause := (WT_D2_1A_pause + WT_D2_2A_pause) / 2]
  sc_dt_UTR[, WT_D4_pause := (WT_D4_1A_pause + WT_D4_2A_pause) / 2]
  plot <- ggplot(data = sc_dt_UTR[WT_D0_1A_rpc >= 1 & WT_D0_2A_rpc >= 1 &
                         WT_D2_1A_rpc >= 1 & WT_D2_2A_rpc >= 1 &
                         WT_D4_1A_rpc >= 1 & WT_D4_2A_rpc >= 1 &
                         WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                         WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
                         WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & length >= 300]) + xlim(-50, 200) +
    stat_summary(aes(position, WT_D0_pause, color = '#1F78B4'), fun.y = "mean", geom = "line", size= 1.25) +
    stat_summary(aes(position, WT_D2_pause, color = '#999999'), fun.y = "mean", geom = "line", size= 1.25) +
    stat_summary(aes(position, WT_D4_pause, color = '#E31A1C'), fun.y = "mean", geom = "line", size= 1.25) +
    scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("#1F78B4", "#999999", "#E31A1C"), name = "")
  plot <- plot + theme_classic(20) + labs(y = "Norm. ribosome occupancy", x = "Codon position") +
    theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
          panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
  ggsave("/Users/KevinStein/Desktop/startMetagene_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

  ### Metagene from stop
  sc_dt_UTR[, length := (length(position) - 133-66-1), by = orf] # subtract flanking 7 codons and stop codon
  sc_dt_UTR[, stopdist := position - (length + 1)]

  plot <- ggplot(data = sc_dt_UTR[WT_D0_1A_rpc >= 1 & WT_D0_2A_rpc >= 1 &
                            WT_D2_1A_rpc >= 1 & WT_D2_2A_rpc >= 1 &
                            WT_D4_1A_rpc >= 1 & WT_D4_2A_rpc >= 1 &
                            WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                            WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
                            WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & length >= 300]) + xlim(-200, 50) +
    stat_summary(aes(stopdist, WT_D4_pause, color = '#E31A1C'), fun.y = "mean", geom = "line", size= 1.25) +
    stat_summary(aes(stopdist, WT_D2_pause, color = '#999999'), fun.y = "mean", geom = "line", size= 1.25) +
    stat_summary(aes(stopdist, WT_D0_pause, color = '#1F78B4'), fun.y = "mean", geom = "line", size= 1.25) +
    scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("#1F78B4", "#999999", "#E31A1C"), name = "")
  plot <- plot + theme_classic(20) + labs(y = "Norm. ribosome occupancy", x = "Codon position") +
    theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
          panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
  ggsave("/Users/KevinStein/Desktop/stopMetagene_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

  
  ### Distribution of pause score by residue ###
  temp <- sc_dtA[WT_D0_1A_rpc >= 1 & WT_D0_2A_rpc >= 1 & 
                   WT_D2_1A_rpc >= 1 & WT_D2_2A_rpc >= 1 &
                   WT_D4_1A_rpc >= 1 & WT_D4_2A_rpc >= 1 &
                   WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 &
                   WT_D2_1A_sum >= 64 & WT_D2_2A_sum >= 64 &
                   WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 &
                   position > 20 & stopdist < -20] 
  
  ggplot(temp, aes(residue, WT_D0_pause)) + geom_violin(fill = "gray40") +
    geom_violin(aes(residue, WT_D4_pause), color = 'red', position = "dodge") + coord_cartesian(ylim = c(0,10))
  
  temp[, residue1 := paste0(residue,"1")]
  plot <- ggplot(temp, aes(residue, WT_D0_pause)) +
    geom_boxplot(aes(residue1, WT_D4_pause), fill = '#F2CF8C', color = '#E7A427', width = 0.8) + 
    geom_boxplot(fill = "gray70", color = "gray40", width = 0.8) + scale_y_log10(breaks = c(0.01,1,100), labels = c(0.01,1,100))
  plot <- plot + theme_minimal(20) + labs(y = "Pause score", x = "") +
    theme(panel.grid.major.x = element_blank())
  ggsave("/Users/KevinStein/Desktop/pauseScoreDistribution.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
  
  
### RQC expression ###
expression <- readRDS("doc/expression.rds")

WT_RQCexpression <- as.data.table(read.csv("/Users/KevinStein/Desktop/WT_RQCexpression.csv", header = TRUE))
names(WT_RQCexpression)    

ggplot(WT_RQCexpression, aes(orf, log2FoldChange)) + geom_col()

RQCexpression <- as.data.table(read.csv("/Users/KevinStein/Desktop/yeast_RQCexpression.csv", header = TRUE))

plot <- ggplot(RQCexpression, aes(gene, log2FoldChange, fill = ID)) + geom_col(position = "dodge") +
  labs(y = "Fold change (Day4 / Day0), log2", x = "", fill = "") + 
  scale_fill_discrete(limits = c("a", "b"), labels = c("WT", "sch9"))
plot <- plot + theme_classic(15) +
  theme(legend.position = c(.8,.3), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=12),
        axis.line.x = element_blank(), axis.text = element_text(size = 15, color = "black"), axis.ticks.x = element_blank())
ggsave("/Users/KevinStein/Desktop/RQCexpression_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

wormRQCexpression <- as.data.table(read.csv("/Users/KevinStein/Desktop/N2_RQCexpression.csv", header = TRUE))
names(wormRQCexpression)

plot <- ggplot(wormRQCexpression, aes(orf, log2FoldChange)) + geom_col() +
  labs(y = "Fold change (Day12 / Day1), log2", x = "")
plot <- plot + theme_classic(15) +
  theme(axis.line.x = element_blank(), axis.text.y = element_text(size = 12, color = "black"), axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 12, color = 'black'))
ggsave("/Users/KevinStein/Desktop/RQCexpression_worm.pdf", plot, width = 4, height = 6, dpi = 300, useDingbats = F)


### Survival curve ###
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

survival <- as.data.table(read.csv("/Users/KevinStein/Desktop/survival.csv", header = TRUE))
names(survival)
survival_dt <- as.data.table(summarySE(survival, measurevar="per_survival", groupvars=c("day")))
spline_int <- as.data.frame(spline(survival_dt$day, survival_dt$per_survival))

plot <- ggplot(survival_dt, aes(x = day, y = per_survival)) + 
  geom_smooth(color = "gray50", se = FALSE, size = 1.25) + geom_point(color = "black", size = 2.5) +
  geom_errorbar(aes(ymin=per_survival-se, ymax=per_survival+se), width=.2, size = 0.6) +
  labs(y = "Percent survival", x = "Age (days)") + scale_x_continuous(breaks = c(2,4,6,8,10,12,14,16))
plot <- plot + theme_minimal(20) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 15, color = "black"))
ggsave("/Users/KevinStein/Desktop/survivalCurve.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Correlation of datasets ###

install.packages("corrplot")
library(corrplot)
library(tidyr)
names(expression)
M <- cor(expression[, c(35,37,39,41,43,45,23,25,27,29,31,33)])
M <- cor(expression[, c(35,37,39,41,43,45)])
head(round(M, 2))
corrplot(M, method="color",  
         type="upper", 
         tl.col="black", tl.srt=45)

pdf(file = "/Users/KevinStein/Desktop/yeast_correlation.pdf", width = 6.5, height = 5.5, useDingbats = F)
corrplot(M, method="color",  
         type="upper", 
         tl.col="black", tl.srt=45, cl.align.text = 'l', addCoef.col = "white", number.cex = 0.7
         )
dev.off()


### Pausing correlation ###
D0_1 <- NULL
D0_2 <- NULL
D2_1 <- NULL
D2_2 <- NULL
D4_1 <- NULL
D4_2 <- NULL
sch9D0_1 <- NULL
sch9D0_2 <- NULL
sch9D2_1 <- NULL
sch9D2_2 <- NULL
sch9D4_1 <- NULL
sch9D4_2 <- NULL
for (i in aa) {
  D0_1a <- mean(residueMean[residue == i]$WT_D0_1A_pause)
  D0_1 <- c(D0_1, D0_1a)
  D0_2a <- mean(residueMean[residue == i]$WT_D0_2A_pause)
  D0_2 <- c(D0_2, D0_2a)
  D2_1a <- mean(residueMean[residue == i]$WT_D2_1A_pause)
  D2_1 <- c(D2_1, D2_1a)
  D2_2a <- mean(residueMean[residue == i]$WT_D2_2A_pause)
  D2_2 <- c(D2_2, D2_2a)
  D4_1a <- mean(residueMean[residue == i]$WT_D4_1A_pause)
  D4_1 <- c(D4_1, D4_1a)
  D4_2a <- mean(residueMean[residue == i]$WT_D4_2A_pause)
  D4_2 <- c(D4_2, D4_2a)
  sch9D0_1a <- mean(residueMean[residue == i]$sch9_D0_1A_pause)
  sch9D0_1 <- c(sch9D0_1, sch9D0_1a)
  sch9D0_2a <- mean(residueMean[residue == i]$sch9_D0_2A_pause)
  sch9D0_2 <- c(sch9D0_2, sch9D0_2a)
  sch9D2_1a <- mean(residueMean[residue == i]$sch9_D2_1A_pause)
  sch9D2_1 <- c(sch9D2_1, sch9D2_1a)
  sch9D2_2a <- mean(residueMean[residue == i]$sch9_D2_2A_pause)
  sch9D2_2 <- c(sch9D2_2, sch9D2_2a)
  sch9D4_1a <- mean(residueMean[residue == i]$sch9_D4_1A_pause)
  sch9D4_1 <- c(sch9D4_1, sch9D4_1a)
  sch9D4_2a <- mean(residueMean[residue == i]$sch9_D4_2A_pause)
  sch9D4_2 <- c(sch9D4_2, sch9D4_2a)
}
pauseMean <- data.table(aa, D0_1, D0_2, D2_1, D2_2, D4_1, D4_2,
                        sch9D0_1, sch9D0_2, sch9D2_1, sch9D2_2, sch9D4_1, sch9D4_2)

M <- cor(pauseMean[, c(2:7)])
head(round(M, 2))
corrplot(M, method="color",  
         type="lower", 
         tl.col="black")

pdf(file = "/Users/KevinStein/Desktop/yeast_pauseScoreCorrelation.pdf", width = 6.5, height = 5.5, useDingbats = F)
corrplot(M, method="color",  
         type="lower", diag = FALSE,
         tl.col="black", tl.srt=45, cl.align.text = 'l', addCoef.col = "white", number.cex = 1
)
dev.off()


### Polysome ###
polysome <- as.data.table(polysome)
names(polysome)
polysome[, WT_D0_2_norm := (WT_D0_2+1) / (max(na.omit(polysome$WT_D0_2))+1)]
polysome[, WT_D2_2_norm := (WT_D2_2+1) / (max(na.omit(polysome$WT_D2_2))+1)]
polysome[, WT_D4_2_norm := (WT_D4_2+1) / (max(na.omit(polysome$WT_D4_2))+1)]
polysome[, WT_D2_1_norm := (WT_D2_1+1) / (max(na.omit(polysome$WT_D2_1))+1)]
polysome[, WT_D4_1_norm := (WT_D4_1+1) / (max(na.omit(polysome$WT_D4_1))+1)]
polysome[, WT_D0_2l_norm := (WT_D0_2long+1) / (max(na.omit(polysome$WT_D0_2long))+1)]
polysome[, WT_D2_2l_norm := (WT_D2_2long+1) / (max(na.omit(polysome$WT_D2_2long))+1)]
polysome[, WT_D4_2l_norm := (WT_D4_2long+1) / (max(na.omit(polysome$WT_D4_2long))+1)]
polysome[, WT_D0_1_norm := (WT_D0_1+1) / (max(na.omit(polysome$WT_D0_1))+1)]

ggplot(polysome, aes(Position, (WT_D0_2_norm))) + geom_line(size = 1.25, color = "gray30") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/yeastD0_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

ggplot(polysome, aes(Position, (WT_D2_1_norm))) + geom_line(size = 1.25, color = "#1F78B4") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
ggplot(polysome, aes(Position, (WT_D2_2l_norm))) + geom_line(size = 1.25, color = "#1F78B4") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())

ggplot(polysome, aes(Position, (WT_D4_1_norm))) + geom_line(size = 1.25, color = "#E31A1C") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
ggplot(polysome, aes(Position, (WT_D4_2l_norm))) + geom_line(size = 1.25, color = "#E31A1C") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())

  
