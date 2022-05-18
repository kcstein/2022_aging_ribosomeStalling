### Extended Data Figure 2 ###


# Fig S2a: Polysome profiles
polysome <- read.csv("/Users/KevinStein/Desktop/Manuscript/Figures/Plots/polysomeProfilesWorms.csv", stringsAsFactors = T, header = T)

plot <- ggplot(polysome, aes(Position, D1_1)) + geom_line(size = 1.25, color = "gray30") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/D1_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(polysome, aes(Position, D12_1)) + geom_line(size = 1.25, color = "#E31A1C") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/D12_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(polysome, aes(Position, (N2_D6_1+1))) + geom_line(size = 1.25, color = "#1F78B4") + 
  scale_y_log10() + scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(20) + labs(y = "", x = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/D6_polysome.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S2b: Metagene from start
plot <- ggplot(data = N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & N2_D6_2A_rpc >= 0.5 &
                               D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                               D1_1A_sum >= 64 & D1_2A_sum >= 64 & N2_D6_2A_sum >= 64 &
                               D12_1A_sum >= 64 & D12_2A_sum >= 64 & length >= 300,
                             .(position, D1 = ((D1_1A_pause + D1_2A_pause) / 2),
                               D12 = ((D12_1A_pause + D12_2A_pause) / 2))]) + xlim(-7, 200) +
  stat_summary(aes(position, D12, color = '#E7298A'), fun.y = "mean", geom = "line", size=1.25) +
  stat_summary(aes(position, D1, color = 'gray40'), fun.y = "mean", geom = "line", size=1.25) +
  scale_color_manual(labels = c("Day 1", "Day 12"), values = c("gray40", "#E7298A"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Norm. ribosome occupancy", x = "Codon position") +
  theme(legend.position = c(.6,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/startMetagene_worm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S2c: Metagene from stop
plot <- ggplot(data = N2_dtA[D1_1A_rpc >= 0.5 & D1_2A_rpc >= 0.5 & N2_D6_2A_rpc >= 0.5 &
                               D12_1A_rpc >= 0.5 & D12_2A_rpc >= 0.5 & 
                               D1_1A_sum >= 64 & D1_2A_sum >= 64 & N2_D6_2A_sum >= 64 &
                               D12_1A_sum >= 64 & D12_2A_sum >= 64 & length >= 300,
                             .(stopdist, D1 = ((D1_1A_pause + D1_2A_pause) / 2),
                               D12 = ((D12_1A_pause + D12_2A_pause) / 2))]) + xlim(-200, 7) +
  stat_summary(aes(stopdist, D1, color = 'gray40'), fun.y = "mean", geom = "line", size=1.25) +
  stat_summary(aes(stopdist, D12, color = '#E7298A'), fun.y = "mean", geom = "line", size=1.25) +
  scale_color_manual(labels = c("Day 1", "Day 12"), values = c("gray40", "#E7298A"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Norm. ribosome occupancy", x = "Codon position") +
  theme(legend.position = c(.3,.9), legend.justification = c("left", "top"), legend.title=element_blank(), legend.text = element_text(size=16),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/stopMetagene_worm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S2d: Correlation of datasets
library(corrplot)
library(tidyr)
N2_expression_deseq <- readRDS("doc/N2_expression_deseq.rds")
names(N2_expression_deseq)
M <- cor(N2_expression_deseq[, c(18,20,22,24,26,28,30)])
M <- cor(N2_expression_deseq[, c(18,20,30,22,24)])
head(round(M, 2))
corrplot(M, method="color",  
         type="upper", 
         tl.col="black", tl.srt=45)

pdf(file = "/Users/KevinStein/Desktop/worm_correlation.pdf", width = 6.5, height = 5.5, useDingbats = F)
corrplot(M, method="color",  
         type="upper", 
         tl.col="black", tl.srt=45, cl.align.text = 'l', addCoef.col = "white", number.cex = 0.9
)
dev.off()


# Fig S2e: Gene-level expression
N2_expression_deseq <- readRDS("doc/N2_expression_deseq.rds")

plot <- ggplot(N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 0]) + 
  geom_point(aes(log2FoldChange, -log10(padj)), fill = "gray75", color = "gray75", alpha = 0.5, size = 2) + 
  geom_point(data = N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 1], aes(log2FoldChange, -log10(padj)), color = factor(group), fill = "#E31A1C", alpha = 0.5, size = 2) + 
  geom_point(data = N2_expression_deseq[D1_1A_sum > 64 & D1_2A_sum > 64 & D12_1A_sum > 64 & D12_2A_sum > 64 & group == 2], aes(log2FoldChange, -log10(padj), color = factor(group)), fill = "#377EB8", alpha = 0.5, size = 2) + 
  scale_color_manual(limits = c("1","2"), labels = c(">2 fold up",">2 fold down"), values = c("#9970AB", "#5AAE61"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "log2 (fold change)") +
  theme(legend.position = c(.01,1.1), legend.justification = c("left", "top"), legend.title=element_text(size=16), legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Deseq.gene.expression_volcano_worms.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# FigS2f: Gene Ontology
GO <- read.csv("/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/RiboSeq_Analysis/GO_translatome.csv", stringsAsFactors = T, header = T)
GO <- as.data.table(GO)
plot <- ggplot(GO[Organism == "worm"], aes(x = reorder(Description, -foldchange), y = foldchange, fill = factor(Group))) + geom_col(position = "dodge", color = "black", size = 0.2) + 
  coord_flip() +
  scale_fill_manual(limits = c("down","up"), values = c("#5AAE61", "#9970AB"), name = "")
plot <- plot + theme_classic(12) + labs(y = "fold enrichment", x = "") +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text = element_text(color = "black", size = 10),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/GO_BP_translatomeWorm.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S2g: Ribosomal protein rps-3
View(N2_dtA[position == 1 & WBgene %in% c("WBGene00010624", "WBGene00017799", "WBGene00020348", "WBGene00004494", "WBGene00004472", "WBGene00004471", "WBGene00001340", "WBGene00004326", "WBGene00004431",
                          "WBGene00016740", "WBGene00004453", "WBGene00007859", "WBGene00004490", "WBGene00004454", "WBGene00004452", "WBGene00004425", "WBGene00004410", "WBGene00004355",
                          "WBGene00017816", "WBGene00006794", "WBGene00017245", "WBGene00010097", "WBGene00004469", "WBGene00000479", "WBGene00004482", "WBGene00000474", "WBGene00004441",
                          "WBGene00015133", "WBGene00012556", "WBGene00004434", "WBGene00015487", "WBGene00016249", "WBGene00004415", "WBGene00004433", "WBGene00004918", "WBGene00004435",
                          "WBGene00004917", "WBGene00015185", "WBGene00004916", "WBGene00012484", "WBGene00016142", "WBGene00004413", "WBGene00015461", "WBGene00016653", "WBGene00016493", "WBGene00022739", "WBGene00004419"), c(70,71)])

plot <- ggplot(data=N2_dtA[WBgene == "WBGene00004472", .(position, D1 = movingAverage(((D1_1A_rpm + D1_2A_rpm) / 2), n=5, center=T), 
                                                 D6 = movingAverage(N2_D6_2A_rpm, n=5, center=T),
                                               D12 = movingAverage(((D12_1A_rpm + D12_2A_rpm) / 2), n=5, center=T))]) + 
  geom_line(aes(position, D1), color = "gray40", size = 1.25) +
  geom_line(aes(position, D6), color = "#F28CC1", size = 1.25) +
  geom_line(aes(position, D12), color = "#E7298A", size = 1.25) +
  scale_x_continuous(expand = expand_scale())
plot <- plot + theme_classic(16) + labs(y = "Pause score", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 16, color = "black"))
ggsave("/Users/KevinStein/Desktop/rps-3.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)



# Fig S2h: Residue pausing
plot <- ggplot(pauseMean, aes(D1, D12)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_point(size = 3, color = "#E7298A", alpha = 0.6) +
  xlim(0.85,1.2) + ylim(0.85,1.2)
plot <- plot + theme_classic(20) + labs(y = "Day 12 pause score", x = "Day 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_wormsPoints.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Fig S2i: Residue pausing by replicates
plot <- ggplot(pauseMean, aes(D1_1, D1_2, color = aa)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) +
  xlim(0.85,1.25) + ylim(0.85,1.25)
plot <- plot + theme_classic(20) + labs(y = "Day 1, rep 2 pause score", x = "Day 1, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day1Replicates_wormsResidues.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(pauseMean, aes(D12_1, D12_2, color = aa)) +  
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed", size = 1) +
  geom_text(label = aa, size = 5) +
  xlim(0.85,1.25) + ylim(0.85,1.25)
plot <- plot + theme_classic(20) + labs(y = "Day12, rep2 pause score", x = "Day 12, rep 1 pause score") +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/allSites_Day12Replicates_wormsResidues.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

