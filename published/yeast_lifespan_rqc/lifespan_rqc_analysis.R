### Lifespan and stalling correlation ###
stall <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Chronological/Expt05_Lifespan/Data/GlobalAnalysis_StallingReporter/Analysis/Brandman.csv", header = T, stringsAsFactors = T))
lifespan <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Chronological/Expt05_Lifespan/Data/GlobalAnalysis_StallingReporter/Analysis/Kaeberlein.csv", header = T, stringsAsFactors = T))
setkeyv(stall, c("ORF"))
setkeyv(lifespan, c("ORF"))


lifespan[, summed_Surv := X2.week + X5.week + X7.weeks]
combined <- lifespan[stall]
combined <- combined[!is.na(GFP_R12) & !is.na(summed_Surv)]
quantile(combined$Norm, c(0.75), na.rm = TRUE)
summary(combined$GFP_R12)
min(na.omit(combined$summed_Surv))
combined[, Norm_adj := summed_Surv + 1.24] # make all positive values
summary(combined$Norm_adj)
combined[, foldchange := 2^GFP_R12]
combined[, lifespan := log2(Norm_adj)]


mean(stall$foldchange)
sd(stall$foldchange)
mean(lifespan$Norm)
sd(lifespan$Norm)


plot <- ggplot(combined[foldchange < x+y], aes("<1sd", log10(Norm_adj))) + 
  geom_hline(yintercept = log10(mean(combined$Norm_adj)), size = 1, color = "gray75", linetype = "dashed") +
  geom_violin(fill = "#E5F5E0") + geom_boxplot(fill = "#E5F5E0", width = 0.3) +
  geom_violin(data = combined[(foldchange > x+y) & (foldchange < x+3*y)], aes("1sd", log10(Norm_adj)), fill = "#A1D99B") +
  geom_boxplot(data = combined[(foldchange > x+y) & (foldchange < x+3*y)], aes("1sd", log10(Norm_adj)), fill = "#A1D99B", width = 0.3) +
  geom_violin(data = combined[!is.na(Norm) & foldchange > x+3*y], aes("3sd", log10(Norm_adj)), fill = "#41AB5D") +
  geom_boxplot(data = combined[!is.na(Norm) & foldchange > x+3*y], aes("3sd", log10(Norm_adj)), fill = "#41AB5D", width = 0.3) +
  geom_point(data = combined[ORF == "YMR247C" | ORF == "YPL009C" | ORF == "YDR333C" | ORF == "YDR266C" | ORF == "YMR116C"], aes("3sd", log10(Norm_adj)), color = '#E7298A', size = 2) +
  geom_text(data = combined[ORF == "YPL009C" | ORF == "YDR266C"], aes("3sd", log10(Norm_adj), label = ORF), color = '#E7298A', nudge_x = 0.25, nudge_y = 0.08) +
  geom_text(data = combined[ORF == "YDR333C"], aes("3sd", log10(Norm_adj), label = ORF), color = '#E7298A', nudge_x = 0.25, nudge_y = 0.08) +
  geom_text(data = combined[ORF == "YMR247C"], aes("3sd", log10(Norm_adj), label = ORF), color = '#E7298A', nudge_x = 0.25, nudge_y = -0.085)
plot <- plot + theme_classic(20) + labs(y = "Norm. Lifespan", x = "GFP-R12 abundance") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Lifespan_Stalling_all.pdf", plot, width = 4, height = 6, dpi = 300, useDingbats = F)



# Examine RQC flux of various groups of strains
tor <- read.csv("/Users/KevinStein/Desktop/Manuscript/Data/Resubmission Data/TOR_related_genes.csv", header = T, stringsAsFactors = T)
RLS <- read.csv("/Users/KevinStein/Desktop/Manuscript/Data/datasets/RA_short.csv", header = T, stringsAsFactors = T)
tor <- as.data.table(tor)
RLS <- as.data.table(RLS)
setkeyv(tor, c("gene"))
setkeyv(RLS, c("gene"))

x1 <- mean(combined$Norm_adj)
y1 <- sd(combined$Norm_adj)
x1+y1

plot <- ggplot(data = combined[Norm_adj > x1+y1], aes("2CLS", log2(foldchange))) + geom_boxplot(fill = "gray60", notch = T) +
  geom_boxplot(data = combined[foldchange > x+3*y], aes("1RQC", log2(foldchange)), fill = "#B2DF8A", notch = T) +
  geom_boxplot(data = stall[stall$ORF %in% tor$gene], aes("3TOR", log2(foldchange)), fill = "gray60", notch = T) +
  geom_boxplot(data = stall[stall$ORF %in% RLS$gene], aes("5RLS", log2(foldchange)), fill = "gray60", notch = T) +
  geom_boxplot(data = combined[(!combined$ORF %in% tor$gene) & Norm_adj > x1+y1], aes("4CLS_noTOR", log2(foldchange)), fill = "gray60", notch = T) +
  geom_point(data = stall[ORF == "YHR205W"], aes("2CLS", log2(foldchange)), color = '#E7298A', size = 2) +
  labs(y = "RQC flux (GFP-R12 fold change, log2)", x = "")
plot <- plot + theme_classic(15) + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/RQCflux.pdf", plot, width = 5, height = 5, dpi = 300, useDingbats = F)

length(combined[Norm_adj > x1+y1]$ORF) #CLS
length(combined[foldchange > x+3*y]$ORF) #RQC
length(stall[stall$ORF %in% tor$gene]$ORF) #TOR
length(stall[stall$ORF %in% RLS$gene]$ORF) #RLS
length(combined[(!combined$ORF %in% tor$gene) & Norm_adj > x1+y1]$ORF) #CLS no TOR
t <- wilcox.test(combined[foldchange > x+3*y]$foldchange, combined[Norm_adj > x1+y1]$foldchange)
t$p.value
t <- wilcox.test(combined[foldchange > x+3*y]$foldchange, stall[stall$ORF %in% tor$gene]$foldchange)
t$p.value
t <- wilcox.test(combined[foldchange > x+3*y]$foldchange, stall[stall$ORF %in% RLS$gene]$foldchange)
t$p.value
t <- wilcox.test(combined[foldchange > x+3*y]$foldchange, combined[(!combined$ORF %in% tor$gene) & Norm_adj > x1+y1]$foldchange)
t$p.value
