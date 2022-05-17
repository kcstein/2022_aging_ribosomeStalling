### Disome motifs ###

# Inhibitory codons
library(data.table)
dicodon <- readRDS("dicodon.rds")
names(dicodon)[names(dicodon) == "position"] <- "codon_position"
names(dicodon)[names(dicodon) == "ID"] <- "name"
dicodon1 <- dicodon[, c(1:2,6)]
dicodon_dt <- sc_dtA[dicodon1, allow.cartesian = TRUE]
dicodon_dt[, codon_position := codon_position + 1]
dicodon_dt[, adjusted := position - codon_position]


### Normalize across polybasic interval ###
test <- dicodon_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("name"))
setkeyv(dicodon_dt, c("name"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = name]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = name]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = name]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = name]
test[, sch9_D0_1A_rpc_adjusted := mean(sch9_D0_1A), by = name]
test[, sch9_D0_2A_rpc_adjusted := mean(sch9_D0_2A), by = name]
test[, sch9_D4_1A_rpc_adjusted := mean(sch9_D4_1A), by = name]
test[, sch9_D4_2A_rpc_adjusted := mean(sch9_D4_2A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(dicodon_dt$name, test$name))
dicodon_dt <- cbind(dicodon_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
dicodon_dt <- cbind(dicodon_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
dicodon_dt <- cbind(dicodon_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
dicodon_dt <- cbind(dicodon_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
dicodon_dt <- cbind(dicodon_dt, sch9_D0_1A_rpc_adjusted = test[i]$sch9_D0_1A_rpc_adjusted)
dicodon_dt <- cbind(dicodon_dt, sch9_D0_2A_rpc_adjusted = test[i]$sch9_D0_2A_rpc_adjusted)
dicodon_dt <- cbind(dicodon_dt, sch9_D4_1A_rpc_adjusted = test[i]$sch9_D4_1A_rpc_adjusted)
dicodon_dt <- cbind(dicodon_dt, sch9_D4_2A_rpc_adjusted = test[i]$sch9_D4_2A_rpc_adjusted)
dicodon_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
dicodon_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
dicodon_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
dicodon_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
dicodon_dt[, sch9_D0_1A_norm := sch9_D0_1A / sch9_D0_1A_rpc_adjusted]
dicodon_dt[, sch9_D0_2A_norm := sch9_D0_2A / sch9_D0_2A_rpc_adjusted]
dicodon_dt[, sch9_D4_1A_norm := sch9_D4_1A / sch9_D4_1A_rpc_adjusted]
dicodon_dt[, sch9_D4_2A_norm := sch9_D4_2A / sch9_D4_2A_rpc_adjusted]
#saveRDS(dicodon_dt, "dicodon_dt.rds")

dicodon_dt <- readRDS("dicodon_dt.rds")
plot <- ggplot(data = dicodon_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                                  WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                                  codon_position > 25 & (codon_position - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/dicodon.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = dicodon_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                           WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                           codon_position > 25 & (codon_position - length) < -25,
                         .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                           WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.3))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/dicodon_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(data = dicodon_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                                  WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                                  sch9_D0_1A_rpc_adjusted >= 0.1 & sch9_D0_2A_rpc_adjusted >= 0.1 & 
                                  sch9_D4_1A_rpc_adjusted >= 0.1 & sch9_D4_2A_rpc_adjusted >= 0.1 & 
                                  codon_position > 25 & (codon_position - length) < -25,
                                .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                                  WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T),
                                  S0 = movingAverage(((sch9_D0_1A_norm + sch9_D0_2A_norm) / 2), n=1, center=T),
                                  S4 = movingAverage(((sch9_D4_1A_norm + sch9_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, S0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray60') +
  stat_summary(aes(adjusted, S4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#00CED1') +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, S0), fun.y = "mean", geom = "line", size=1.25, color = 'gray60') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') + 
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, S4), fun.y = "mean", geom = "line", size=1.25, color = '#00CED1') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/dicodon_sch9.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Examples
plot <- ggplot(data=sc_dtA[orf == "YGR186W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_rect(xmin = 298, xmax = 338, ymin = -Inf, ymax = Inf, fill = "#F2CF8C", alpha = 0.1) +
  geom_vline(xintercept = 318, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YGR186W_GGGaa318.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YGR186W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 318, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(298,338))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/YGR186W_GGGaa318_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(data=sc_dtA[orf == "YPR133C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_rect(xmin = 211, xmax = 251, ymin = -Inf, ymax = Inf, fill = "#F2CF8C", alpha = 0.1) +
  geom_vline(xintercept = 231, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/SPN1_PDGaa231.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YPR133C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 231, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(211,251))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/SPN1_PDGaa231_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YPL127C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_rect(xmin = 143, xmax = 183, ymin = -Inf, ymax = Inf, fill = "#F2CF8C", alpha = 0.1) +
  geom_vline(xintercept = 163, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/HHO1_KKKaa163.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YPL127C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 163, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(143,183))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/HHO1_KKKaa163_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YLR196W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_rect(xmin = 257, xmax = 297, ymin = -Inf, ymax = Inf, fill = "#F2CF8C", alpha = 0.1) +
  geom_vline(xintercept = 277, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/PWP1_KKKKKKSKaa277.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YLR196W", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 277, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(257,297))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/PWP1_KKKKKKSKaa277_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YPR093C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_rect(xmin = 207, xmax = 247, ymin = -Inf, ymax = Inf, fill = "#F2CF8C", alpha = 0.1) +
  geom_vline(xintercept = 227, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/ASR1_RRaa227.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data=sc_dtA[orf == "YPR093C", .(position, WT0 = movingAverage(((WT_D0_1A_pause + WT_D0_2A_pause) / 2), n=4, center=T), 
                                               WT4 = movingAverage(((WT_D4_1A_pause + WT_D4_2A_pause) / 2), n=4, center=T))]) + 
  geom_vline(xintercept = 227, color = "gray75", linetype = "longdash", size = 1) +
  geom_line(aes(position, WT4), color = "#E7A427", size = 1.25) +
  geom_line(aes(position, WT0), color = "gray40", size = 1.25) +
  scale_x_continuous(expand = expand_scale(), limits = c(207,247))
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/ASR1_RRaa227_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Other motifs
RKK <- sc_dtA[motif3 == "RKK"]
#RKK <- readRDS("RKK.rds")
names(RKK)[names(RKK) == "position"] <- "codon_position"
names(RKK)[names(RKK) == "ID"] <- "name"
RKK1 <- RKK[, c(1:2,6)]
RKK_dt <- sc_dtA[RKK1, allow.cartesian = TRUE]
RKK_dt[, adjusted := position - codon_position]


### Normalize across polybasic interval ###
test <- RKK_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("name"))
setkeyv(RKK_dt, c("name"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = name]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = name]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = name]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = name]
test[, sch9_D0_1A_rpc_adjusted := mean(sch9_D0_1A), by = name]
test[, sch9_D0_2A_rpc_adjusted := mean(sch9_D0_2A), by = name]
test[, sch9_D4_1A_rpc_adjusted := mean(sch9_D4_1A), by = name]
test[, sch9_D4_2A_rpc_adjusted := mean(sch9_D4_2A), by = name]
test <- test[, .SD[which.min(position)], by = name]
i <- cbind(match(RKK_dt$name, test$name))
RKK_dt <- cbind(RKK_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
RKK_dt <- cbind(RKK_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
RKK_dt <- cbind(RKK_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
RKK_dt <- cbind(RKK_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
RKK_dt <- cbind(RKK_dt, sch9_D0_1A_rpc_adjusted = test[i]$sch9_D0_1A_rpc_adjusted)
RKK_dt <- cbind(RKK_dt, sch9_D0_2A_rpc_adjusted = test[i]$sch9_D0_2A_rpc_adjusted)
RKK_dt <- cbind(RKK_dt, sch9_D4_1A_rpc_adjusted = test[i]$sch9_D4_1A_rpc_adjusted)
RKK_dt <- cbind(RKK_dt, sch9_D4_2A_rpc_adjusted = test[i]$sch9_D4_2A_rpc_adjusted)
RKK_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
RKK_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
RKK_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
RKK_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
RKK_dt[, sch9_D0_1A_norm := sch9_D0_1A / sch9_D0_1A_rpc_adjusted]
RKK_dt[, sch9_D0_2A_norm := sch9_D0_2A / sch9_D0_2A_rpc_adjusted]
RKK_dt[, sch9_D4_1A_norm := sch9_D4_1A / sch9_D4_1A_rpc_adjusted]
RKK_dt[, sch9_D4_2A_norm := sch9_D4_2A / sch9_D4_2A_rpc_adjusted]
#saveRDS(RKK_dt, "RKK_dt.rds")

RKK_dt <- readRDS("RKK_dt.rds")
plot <- ggplot(data = RKK_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/RKK.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = RKK_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.2))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/RKK_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(data = RKK_dt[WT_D0_1A_rpc_adjusted >= 0.1 & WT_D0_2A_rpc_adjusted >= 0.1 & 
                               WT_D4_1A_rpc_adjusted >= 0.1 & WT_D4_2A_rpc_adjusted >= 0.1 & 
                               sch9_D0_1A_rpc_adjusted >= 0.1 & sch9_D0_2A_rpc_adjusted >= 0.1 & 
                               sch9_D4_1A_rpc_adjusted >= 0.1 & sch9_D4_2A_rpc_adjusted >= 0.1 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T),
                               S0 = movingAverage(((sch9_D0_1A_norm + sch9_D0_2A_norm) / 2), n=1, center=T),
                               S4 = movingAverage(((sch9_D4_1A_norm + sch9_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, S0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray60') +
  stat_summary(aes(adjusted, S4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#00CED1') +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, S0), fun.y = "mean", geom = "line", size=1.25, color = 'gray60') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') + 
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, S4), fun.y = "mean", geom = "line", size=1.25, color = '#00CED1') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/RKK_sch9.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


PPG_dt <- readRDS("PPG_dt.rds")
plot <- ggplot(data = PPG_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/PPG.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = PPG_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.2))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/PPG_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


GGG_dt <- readRDS("GGG_dt.rds")
plot <- ggplot(data = GGG_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/GGG.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = GGG_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                               WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 & 
                               codon_position > 25 & (codon_position - length) < -25,
                             .(adjusted, WT0 = movingAverage(((WT_D0_1A_norm + WT_D0_2A_norm) / 2), n=1, center=T),
                               WT4 = movingAverage(((WT_D4_1A_norm + WT_D4_2A_norm) / 2), n=1, center=T))]) + 
  geom_vline(xintercept = 0, size = 1, color = "gray75", linetype = "dashed") +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.2))
plot <- plot + theme_classic(16) + labs(y = "Norm. ribosome occupancy", x = "Distance from inhibitory codon pair (codons)") +
  theme(legend.position = c(.01,.99), legend.title = element_blank(), legend.justification = c("left", "top"), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/GGG_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)




### Metagene around disome positions ###
disome_positions <- as.data.table(read.csv("/Users/KevinStein/Desktop/Manuscript/Data/Resubmission Data/Collisions/disome_positions.csv"))
disome_peaks_dt <- sc_dtA[disome_positions, allow.cartesian = TRUE]
disome_peaks_dt[, adjusted := position - disome]

test <- disome_peaks_dt[adjusted >= -40 & adjusted <= 40]
setkeyv(test, c("disome_ID"))
setkeyv(disome_peaks_dt, c("disome_ID"))
test[, WT_D0_1A_rpc_adjusted := mean(WT_D0_1A), by = disome_ID]
test[, WT_D0_2A_rpc_adjusted := mean(WT_D0_2A), by = disome_ID]
test[, WT_D2_1A_rpc_adjusted := mean(WT_D2_1A), by = disome_ID]
test[, WT_D2_2A_rpc_adjusted := mean(WT_D2_2A), by = disome_ID]
test[, WT_D4_1A_rpc_adjusted := mean(WT_D4_1A), by = disome_ID]
test[, WT_D4_2A_rpc_adjusted := mean(WT_D4_2A), by = disome_ID]
test <- test[, .SD[which.min(position)], by = disome_ID]
i <- cbind(match(disome_peaks_dt$disome_ID, test$disome_ID))
disome_peaks_dt <- cbind(disome_peaks_dt, WT_D0_1A_rpc_adjusted = test[i]$WT_D0_1A_rpc_adjusted)
disome_peaks_dt <- cbind(disome_peaks_dt, WT_D0_2A_rpc_adjusted = test[i]$WT_D0_2A_rpc_adjusted)
disome_peaks_dt <- cbind(disome_peaks_dt, WT_D2_1A_rpc_adjusted = test[i]$WT_D2_1A_rpc_adjusted)
disome_peaks_dt <- cbind(disome_peaks_dt, WT_D2_2A_rpc_adjusted = test[i]$WT_D2_2A_rpc_adjusted)
disome_peaks_dt <- cbind(disome_peaks_dt, WT_D4_1A_rpc_adjusted = test[i]$WT_D4_1A_rpc_adjusted)
disome_peaks_dt <- cbind(disome_peaks_dt, WT_D4_2A_rpc_adjusted = test[i]$WT_D4_2A_rpc_adjusted)
disome_peaks_dt[, WT_D0_1A_norm := WT_D0_1A / WT_D0_1A_rpc_adjusted]
disome_peaks_dt[, WT_D0_2A_norm := WT_D0_2A / WT_D0_2A_rpc_adjusted]
disome_peaks_dt[, WT_D2_1A_norm := WT_D2_1A / WT_D2_1A_rpc_adjusted]
disome_peaks_dt[, WT_D2_2A_norm := WT_D2_2A / WT_D2_2A_rpc_adjusted]
disome_peaks_dt[, WT_D4_1A_norm := WT_D4_1A / WT_D4_1A_rpc_adjusted]
disome_peaks_dt[, WT_D4_2A_norm := WT_D4_2A / WT_D4_2A_rpc_adjusted]
disome_peaks_dt[, WT_D0_norm := (WT_D0_1A_norm + WT_D0_2A_norm) / 2]
disome_peaks_dt[, WT_D2_norm := (WT_D2_1A_norm + WT_D2_2A_norm) / 2]
disome_peaks_dt[, WT_D4_norm := (WT_D4_1A_norm + WT_D4_2A_norm) / 2]
# saveRDS(disome_peaks_dt, "disome_peaks_dt.rds")

disome_peaks_dt <- readRDS("disome_peaks_dt.rds")
plot <- ggplot(data = disome_peaks_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                                #WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                                WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 &
                                disome > 25 & (disome - length) < -25,
                              .(adjusted, WT0 = movingAverage(WT_D0_norm, n=1, center=T),
                                WT2 = movingAverage(WT_D2_norm, n=1, center=T),
                                WT4 = movingAverage(WT_D4_norm, n=1, center=T))]) +
  # stat_summary(aes(adjusted, WT2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
  #              fun.args=list(conf.int=0.5), fill = '#F2CF8C') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  #stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = '#F2CF8C') +  
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray40", "#F2CF8C", "#E7A427"), name = "") +
  scale_x_continuous(limits = c(-25, 25), expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from disome site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_disomepositions_Age.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = disome_peaks_dt[WT_D0_1A_rpc_adjusted >= 0.5 & WT_D0_2A_rpc_adjusted >= 0.5 & 
                                        #WT_D2_1A_rpc_adjusted >= 1 & WT_D2_2A_rpc_adjusted >= 1 &
                                        WT_D4_1A_rpc_adjusted >= 0.5 & WT_D4_2A_rpc_adjusted >= 0.5 &
                                        disome > 25 & (disome - length) < -25,
                                      .(adjusted, WT0 = movingAverage(WT_D0_norm, n=1, center=T),
                                        WT2 = movingAverage(WT_D2_norm, n=1, center=T),
                                        WT4 = movingAverage(WT_D4_norm, n=1, center=T))]) +
  # stat_summary(aes(adjusted, WT2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
  #              fun.args=list(conf.int=0.5), fill = '#F2CF8C') +
  stat_summary(aes(adjusted, WT4), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#E7A427') +
  stat_summary(aes(adjusted, WT0), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray40') +
  stat_summary(aes(adjusted, WT4), fun.y = "mean", geom = "line", size=1.25, color = '#E7A427') +
  #stat_summary(aes(adjusted, WT2), fun.y = "mean", geom = "line", size=1.25, color = '#F2CF8C') +  
  stat_summary(aes(adjusted, WT0), fun.y = "mean", geom = "line", size=1.25, color = 'gray40') +
  scale_color_manual(labels = c("Day 0","Day 2", "Day 4"), values = c("gray40", "#F2CF8C", "#E7A427"), name = "") +
  scale_x_continuous(limits = c(-15, -5), expand = expand_scale()) + coord_cartesian(ylim = c(0.9,1.3))
plot <- plot + theme_classic(20) + labs(y = "Pause score", x = "Distance from disome site (codons)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stallingMetagene_disomepositions_Age_inset.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
