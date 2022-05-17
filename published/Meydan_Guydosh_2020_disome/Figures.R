# Run MotifAnalysis script

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


### Examples: Ytm1, Hat2, His3
plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YOR272W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=4, center=T),
                                                                            occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=4, center=T))]) + 
  geom_rect(xmin = 242, xmax = 282, ymin = -Inf, ymax = Inf, fill = "#AFEEEE", alpha = 0.01) +
  geom_vline(xintercept = 262, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  geom_line(color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/Ytm1_disomeMonosome.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YOR272W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=3, center=T),
                                                                            occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=3, center=T))]) + 
  geom_vline(xintercept = 262, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  geom_line(color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(242,282)) + coord_cartesian(ylim = c(0,1.5))
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/Ytm1_disomeMonosome_inset.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YEL056W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=4, center=T),
                                                                            occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=4, center=T))]) + 
  geom_rect(xmin = 168, xmax = 208, ymin = -Inf, ymax = Inf, fill = "#AFEEEE", alpha = 0.01) +
  geom_vline(xintercept = 188, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  geom_line(color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/Hat2_disomeMonosome.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YEL056W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=3, center=T),
                                                                            occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=3, center=T))]) + 
  geom_vline(xintercept = 188, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  geom_line(color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(168,208)) + coord_cartesian(ylim = c(0,1.5))
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/Hat2_disomeMonosome_inset.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YOR202W", .(position, occupancy = movingAverage(((WT_mono1_A_rpm + WT_mono2_A_rpm) / 2), n=4, center=T),
                                                                            occupancy2 = movingAverage(((WT_di1_A_rpm + WT_di2_A_rpm + WT_di3_A_rpm) / 3), n=4, center=T))]) + 
  geom_line(color = "gray40", size = 1.25) +
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  scale_x_continuous(expand = expand_scale())
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/His3_disomeMonosome.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(aes(position, occupancy), data=disome_dt[orf == "YOR202W", .(position, occupancy = movingAverage(WT_His3_di_A_rpm, n=4, center=T),
                                                                            occupancy2 = movingAverage(WT_His3_tri_A_rpm, n=4, center=T))]) + 
  geom_line(aes(position, occupancy2), color = '#00CED1', size = 1.25) +
  geom_line(color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
G <- plot + theme_classic(16) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/His3_disomeTrisome.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


