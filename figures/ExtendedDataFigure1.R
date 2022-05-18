### Extended Data Figure 1 ###


# Fig S1a: survival curve
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


# Fig S1b: Metagene from start
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


# Fig S1c: Metagene from stop
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


# Fig S1d: Polysome profiles
polysome <- read.csv("/Users/KevinStein/Desktop/Manuscript/1st Submission/Data/Plots/FigS1/polysomeProfilesYeast.csv", stringsAsFactors = T, header = T)

plot <- ggplot(polysome, aes(Position, (WT_D0_2 + 1))) + geom_line(size = 1.25, color = "gray30") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/yeastD0_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(polysome, aes(Position, (WT_D2_2 + 1))) + geom_line(size = 1.25, color = "#1F78B4") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/yeastD2_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(polysome, aes(Position, (WT_D4_2 + 1))) + geom_line(size = 1.25, color = "#E31A1C") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/yeastD4_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S1e: Correlation of datasets
library(corrplot)
library(tidyr)
names(expression)
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


# Fig S1f: Gene-level expression
plot <- ggplot(WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 0]) + 
  geom_point(aes(log2FoldChange, -log10(padj)), fill = "gray75", color = "gray75", alpha = 0.5, size = 2) + 
  geom_point(data = WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 1], aes(log2FoldChange, -log10(padj), color = factor(group)), fill = "#E31A1C", alpha = 0.5, size = 2) + 
  geom_point(data = WT_expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64 & WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64 & group == 2], aes(log2FoldChange, -log10(padj), color = factor(group)), fill = "#377EB8", alpha = 0.5, size = 2) + 
  scale_color_manual(limits = c("1","2"), labels = c(">2 fold up",">2 fold down"), values = c("#9970AB", "#5AAE61"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "log2 (fold change)") +
  theme(legend.position = c(.01,1.1), legend.justification = c("left", "top"), legend.title=element_text(size=16), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Deseq.gene.expression_volcano_yeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S1g: Gene ontology
GO <- read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/GO_translatome.csv", stringsAsFactors = T, header = T)
GO <- as.data.table(GO)
plot <- ggplot(GO[Organism == "yeast"], aes(x = reorder(Description, -foldchange), y = foldchange, fill = factor(Group))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("down","up"), values = c("#1F78B4", "#E31A1C"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_BP_translatomeYeast.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S1h: Ribosomal protein RPL3
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


# Fig S1i: GCN4
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


# Fig S1j: Residue pausing
plot <- ggplot(pauseMean, aes(D0, D4)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 3, color = "#E7A427", alpha = 0.6) +
  xlim(0.5,1.7) + ylim(0.5,1.7)
plot <- plot + theme_classic(20) + labs(y = "Day 4 pause score", x = "Day 0 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_yeastPoints.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
cor(pauseMean$D0, pauseMean$D4)


# Fig S1k: Residue pausing by replicates
plot <- ggplot(pauseMean, aes(D0_1, D0_2, color = aa)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) +
  xlim(0.5,1.7) + ylim(0.5,1.7)
plot <- plot + theme_classic(20) + labs(y = "Day0, rep2 pause score", x = "Day 0, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day0Replicates_yeastResidues.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMean, aes(D4_1, D4_2, color = aa)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) +
  xlim(0.5,1.7) + ylim(0.5,1.7)
plot <- plot + theme_classic(20) + labs(y = "Day4, rep2 pause score", x = "Day 4, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day4Replicates_yeastResidues.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

