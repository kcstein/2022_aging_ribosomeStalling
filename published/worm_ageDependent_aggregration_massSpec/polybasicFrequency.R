library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

### Grab AA composition ###
ce_comp <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/proteome_AAcomposition_fraction.csv', header = TRUE, stringsAsFactors = TRUE))

dir <- "/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/composition/"
files <- list.files(dir, pattern = "*.csv")
samples <- sub("\\.csv", '', files)

for(i in 1:length(files)) {
  print(files[i])
  j <- as.data.table(lapply(paste0(dir, files[i]), read.csv, header = FALSE, stringsAsFactors = TRUE))
  setnames(j, c('orf', 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'))
  setkeyv(j, c("orf"))
  assign(samples[i], j)
}


### Polybasic analysis ###
dir <- "/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/polybasic/"
files <- list.files(dir, pattern = "*.csv")
samples <- sub("\\.csv", '', files)

for(i in 1:length(files)) {
  print(files[i])
  dt <- as.data.table(read.csv(paste0(dir, files[i]), header = TRUE, stringsAsFactors = TRUE))
  colnames(dt) <- c("orf", "polybasic")
  dt[, name := paste0(orf, "_", polybasic)]
  assign(samples[i], dt)
}

KR6of6 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/KR6of6.csv", header = FALSE, stringsAsFactors = TRUE))
colnames(KR6of6) <- c("orf", "polybasic")
KR6of6[, name := paste0(orf, "_", polybasic)]
KR5of5 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/KR5of5.csv", header = FALSE, stringsAsFactors = TRUE))
colnames(KR5of5) <- c("orf", "polybasic")
KR5of5[, name := paste0(orf, "_", polybasic)]
KR4of4 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/KR4of4.csv", header = FALSE, stringsAsFactors = TRUE))
colnames(KR4of4) <- c("orf", "polybasic")
KR4of4[, name := paste0(orf, "_", polybasic)]
KR6of10 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/KR6of10.csv", header = FALSE, stringsAsFactors = TRUE))
colnames(KR6of10) <- c("orf", "polybasic")
KR6of10[, name := paste0(orf, "_", polybasic)]


length(N2_dtA[position > 0 & stopdist < 0]$orf) #8319550
length(unique(N2_dtA[position > 0 & stopdist < 0]$orf)) #20222
length(N2_dtA[position > 0 & stopdist < 0 & motif2 == "RR"]$orf) #32545
length(unique(N2_dtA[position > 0 & stopdist < 0 & motif2 == "RR"]$orf)) #12264
length(N2_dtA[position > 0 & stopdist < 0 & motif3 == "RRR"]$orf) #3220
length(unique(N2_dtA[position > 0 & stopdist < 0 & motif3 == "RRR"]$orf)) #2226


# Kenyon aggregates: significance of polybasic frequency
k4 <- matrix(c(length(unique(kenyon_aggregates_KR4of4$orf)), length(unique(KR4of4$orf)),
               length(unique(kenyon_aggregates_AAcomposition$orf)) - length(unique(kenyon_aggregates_KR4of4$orf)),
               length(unique(ce_comp$orf)) - length(unique(KR4of4$orf))), nrow = 2)
k4_odds <- fisher.test(k4)
k4_odds$p.value


k2 <- matrix(c(length(unique(kenyon_aggregates_R2of2$orf)), 12264,
               length(unique(kenyon_aggregates_AAcomposition$orf)) - length(unique(kenyon_aggregates_R2of2$orf)),
               20222 - 12264), nrow = 2)
k2_odds <- fisher.test(k2)
k3 <- matrix(c(length(unique(kenyon_aggregates_R3of3$orf)), 2226,
               length(unique(kenyon_aggregates_AAcomposition$orf)) - length(unique(kenyon_aggregates_R3of3$orf)),
               20222 - 2226), nrow = 2)
k3_odds <- fisher.test(k3)


kstall4 <- matrix(c(length(unique(kenyon_aggregatesANDstall_KR4of4$orf)), length(unique(KR4of4$orf)),
               length(unique(kenyon_aggregatesANDstall_AAcomposition$orf)) - length(unique(kenyon_aggregatesANDstall_KR4of4$orf)),
               length(unique(ce_comp$orf)) - length(unique(KR4of4$orf))), nrow = 2)
kstall4_odds <- fisher.test(kstall4)
kstall4_odds$p.value

kNOstall4 <- matrix(c(length(unique(kenyon_aggregatesNOTstall_KR4of4$orf)), length(unique(KR4of4$orf)),
                    length(unique(kenyon_aggregatesNOTstall_AAcomposition$orf)) - length(unique(kenyon_aggregatesNOTstall_KR4of4$orf)),
                    length(unique(ce_comp$orf)) - length(unique(KR4of4$orf))), nrow = 2)
kNOstall4_odds <- fisher.test(kNOstall4)
kNOstall4_odds$p.value


kstall2 <- matrix(c(length(unique(kenyon_aggregatesANDstall_R2of2$orf)), 12264,
               length(unique(kenyon_aggregatesANDstall_AAcomposition$orf)) - length(unique(kenyon_aggregatesANDstall_R2of2$orf)),
               20222 - 12264), nrow = 2)
kstall2_odds <- fisher.test(kstall2)

kNOstall2 <- matrix(c(length(unique(kenyon_aggregatesNOTstall_R2of2$orf)), 12264,
                      length(unique(kenyon_aggregatesNOTstall_AAcomposition$orf)) - length(unique(kenyon_aggregatesNOTstall_R2of2$orf)),
                      20222 - 12264), nrow = 2)
kNOstall2_odds <- fisher.test(kNOstall2)


kstall3 <- matrix(c(length(unique(kenyon_aggregatesANDstall_R3of3$orf)), 2226,
                    length(unique(kenyon_aggregatesANDstall_AAcomposition$orf)) - length(unique(kenyon_aggregatesANDstall_R3of3$orf)),
                    20222 - 2226), nrow = 2)
kstall3_odds <- fisher.test(kstall3)

kNOstall3 <- matrix(c(length(unique(kenyon_aggregatesNOTstall_R3of3$orf)), 2226,
                      length(unique(kenyon_aggregatesNOTstall_AAcomposition$orf)) - length(unique(kenyon_aggregatesNOTstall_R3of3$orf)),
                      20222 - 2226), nrow = 2)
kNOstall3_odds <- fisher.test(kNOstall3)



# Hartl aggregates: significance of polybasic frequency
h4 <- matrix(c(length(unique(hartl_aggregates_KR4of4$orf)), length(unique(KR4of4$orf)),
               length(unique(hartl_aggregates_AAcomposition$orf)) - length(unique(hartl_aggregates_KR4of4$orf)),
               length(unique(ce_comp$orf)) - length(unique(KR4of4$orf))), nrow = 2)
h4_odds <- fisher.test(h4)
h4_odds$p.value

h2 <- matrix(c(length(unique(hartl_aggregates_R2of2$orf)), 12264,
               length(unique(hartl_aggregates_AAcomposition$orf)) - length(unique(hartl_aggregates_R2of2$orf)),
               20222 - 12264), nrow = 2)
h2_odds <- fisher.test(h2)
h3 <- matrix(c(length(unique(hartl_aggregates_R3of3$orf)), 2226,
               length(unique(hartl_aggregates_AAcomposition$orf)) - length(unique(hartl_aggregates_R3of3$orf)),
               20222 - 2226), nrow = 2)
h3_odds <- fisher.test(h3)


hstall4 <- matrix(c(length(unique(hartl_aggregatesANDstall_KR4of4$orf)), length(unique(KR4of4$orf)),
               length(unique(hartl_aggregatesANDstall_AAcomposition$orf)) - length(unique(hartl_aggregatesANDstall_KR4of4$orf)),
               length(unique(ce_comp$orf)) - length(unique(KR4of4$orf))), nrow = 2)
hstall4_odds <- fisher.test(hstall4)
hstall4_odds$p.value

hNOstall4 <- matrix(c(length(unique(hartl_aggregatesNOTstall_KR4of4$orf)), length(unique(KR4of4$orf)),
                      length(unique(hartl_aggregatesNOTstall_AAcomposition$orf)) - length(unique(hartl_aggregatesNOTstall_KR4of4$orf)),
                      length(unique(ce_comp$orf)) - length(unique(KR4of4$orf))), nrow = 2)
hNOstall4_odds <- fisher.test(hNOstall4)
hNOstall4_odds$p.value


hstall2 <- matrix(c(length(unique(hartl_aggregatesANDstall_R2of2$orf)), 12264,
               length(unique(hartl_aggregatesANDstall_AAcomposition$orf)) - length(unique(hartl_aggregatesANDstall_R2of2$orf)),
               20222 - 12264), nrow = 2)
hstall2_odds <- fisher.test(hstall2)

hNOstall2 <- matrix(c(length(unique(hartl_aggregatesNOTstall_R2of2$orf)), 12264,
                    length(unique(hartl_aggregatesNOTstall_AAcomposition$orf)) - length(unique(hartl_aggregatesNOTstall_R2of2$orf)),
                    20222 - 12264), nrow = 2)
hNOstall2_odds <- fisher.test(hNOstall2)


hstall3 <- matrix(c(length(unique(hartl_aggregatesANDstall_R3of3$orf)), 2226,
               length(unique(hartl_aggregatesANDstall_AAcomposition$orf)) - length(unique(hartl_aggregatesANDstall_R3of3$orf)),
               20222 - 2226), nrow = 2)
hstall3_odds <- fisher.test(hstall3)

hNOstall3 <- matrix(c(length(unique(hartl_aggregatesNOTstall_R3of3$orf)), 2226,
                      length(unique(hartl_aggregatesNOTstall_AAcomposition$orf)) - length(unique(hartl_aggregatesNOTstall_R3of3$orf)),
                      20222 - 2226), nrow = 2)
hNOstall3_odds <- fisher.test(hNOstall3)



datasource <- c(rep("Hartl",9), rep("Kenyon",9))
residues <- c("KR4", "KR4", "KR4", "R2", "R2", "R2", "R3", "R3", "R3",
              "KR4", "KR4", "KR4", "R2", "R2", "R2", "R3", "R3", "R3")
stall <- c("a", "n", "y", "a", "n", "y", "a", "n", "y",
           "a", "n", "y", "a", "n", "y", "a", "n", "y")
odds <- c(h4_odds$estimate, hstall4_odds$estimate, hNOstall4_odds$estimate,
          h2_odds$estimate, hstall2_odds$estimate, hNOstall2_odds$estimate,
          h3_odds$estimate, hstall3_odds$estimate, hNOstall3_odds$estimate,
          k4_odds$estimate, kstall4_odds$estimate, kNOstall4_odds$estimate,
          k2_odds$estimate, kstall2_odds$estimate, kNOstall2_odds$estimate,
          k3_odds$estimate, kstall3_odds$estimate, kNOstall3_odds$estimate)
agg_polybasic_dt <- data.table(datasource = datasource,
                               residues = residues,
                               stall = stall,
                               odds = odds)

plot <- ggplot(agg_polybasic_dt[datasource == "Hartl"], aes(x = residues, y = odds, fill = stall)) + 
  geom_col(color = "black", size = 0.5, position = "dodge") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 1) +
  scale_fill_manual(limits = c("a", "n","y"), labels = c("all", "no", "yes"), values = c("gray50", "#FB9A99", "#E31A1C"), name = "Pause site") +
  scale_y_continuous(expand = c(0.0008,0.0008)) +
  scale_x_discrete(labels = c("4 Lys/Arg", "2 Arg", "3 Arg"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio (relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.9,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_HartlAggregation.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)



datasource <- c(rep("Hartl",4), rep("Kenyon",4))
residues <- c("KR4", "KR4", "R2", "R2",
              "KR4", "KR4", "R2", "R2")
stall <- c("a", "y", "a", "y", 
           "a", "y", "a", "y")
odds <- c(h4_odds$estimate, hstall4_odds$estimate, 
          h2_odds$estimate, hstall2_odds$estimate, 
          k4_odds$estimate, kstall4_odds$estimate, 
          k2_odds$estimate, kstall2_odds$estimate)
agg_polybasic_dt <- data.table(datasource = datasource,
                               residues = residues,
                               stall = stall,
                               odds = odds)

plot <- ggplot(agg_polybasic_dt[datasource == "Hartl"], aes(x = residues, y = odds, fill = stall)) + 
  geom_col(color = "black", size = 0.5, position = "dodge") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(limits = c("a","y"), labels = c("all", "yes"), values = c("gray50", "#00CED1"), name = "Pause site") +
  scale_y_continuous(expand = c(0.0008,0.0008)) +
  scale_x_discrete(labels = c("4 Lys/Arg", "2 Arg"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio (relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.1,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_HartlAggregation.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)

plot <- ggplot(agg_polybasic_dt[datasource == "Kenyon"], aes(x = residues, y = odds, fill = stall)) + 
  geom_col(color = "black", size = 0.5, position = "dodge") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(limits = c("a","y"), labels = c("all", "yes"), values = c("gray50", "#00CED1"), name = "Pause site") +
  scale_y_continuous(expand = c(0.0008,0.0008)) +
  scale_x_discrete(labels = c("4 Lys/Arg", "2 Arg"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio (relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.1,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_KenyonAggregation.pdf", plot, width = 4, height = 5, dpi = 300, useDingbats = F)


datasource <- c(rep("Hartl",3), rep("Kenyon",3))
residues <- c("KR4", "R2", "R3", "KR4", "R2", "R3")
odds <- c(h4_odds$estimate, h2_odds$estimate, h3_odds$estimate, 
          k4_odds$estimate, k2_odds$estimate, k3_odds$estimate)
agg_polybasic_dt <- data.table(datasource = datasource,
                               residues = residues,
                               odds = odds)

plot <- ggplot(agg_polybasic_dt[datasource == "Hartl"], aes(x = residues, y = odds)) + 
  geom_col(color = "black", size = 0.5, fill = "#80B1D3") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", size = 1) +
  coord_capped_cart(left = "both") +
  #scale_y_continuous(expand = c(0.0008,0.0008)) +
  scale_x_discrete(labels = c("4 K/R", "2 R", "3 R"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio\n(relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.9,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_HartlAggregation.pdf", plot, width = 3.2, height = 5, dpi = 300, useDingbats = F)

plot <- ggplot(agg_polybasic_dt[datasource == "Kenyon"], aes(x = residues, y = odds)) + 
  geom_col(color = "black", size = 0.5, fill = "#80B1D3") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", size = 1) +
  coord_capped_cart(left = "both", ylim = c(0,1.5)) +
  scale_x_discrete(labels = c("4 K/R", "2 R", "3 R"))
plot <- plot + theme_classic(18) + labs(y = "Odds ratio\n(relative to proteome)", x = "Residue stretch") +
  theme(legend.position = c(0.9,0.9), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/PolybasicEnrichment_KenyonAggregation.pdf", plot, width = 3.2, height = 5, dpi = 300, useDingbats = F)

datasource <- c(rep("Hartl",3), rep("Kenyon",3), rep("Proteome",3))
residues <- c("KR4", "R2", "R3", "KR4", "R2", "R3", "KR4", "R2", "R3")
fraction <- c(length(unique(hartl_aggregates_KR4of4$orf)) / length(unique(hartl_aggregates_AAcomposition$orf)), 
              length(unique(hartl_aggregates_R2of2$orf)) / length(unique(hartl_aggregates_AAcomposition$orf)), 
              length(unique(hartl_aggregates_R3of3$orf)) / length(unique(hartl_aggregates_AAcomposition$orf)), 
              length(unique(kenyon_aggregates_KR4of4$orf)) / length(unique(kenyon_aggregates_AAcomposition$orf)), 
              length(unique(kenyon_aggregates_R2of2$orf)) / length(unique(kenyon_aggregates_AAcomposition$orf)), 
              length(unique(kenyon_aggregates_R3of3$orf)) / length(unique(kenyon_aggregates_AAcomposition$orf)), 
              length(unique(KR4of4$orf)) / length(unique(ce_comp$orf)),
              12264 / 20222,
              2226 / 20222)
agg_polybasic_dt <- data.table(datasource = datasource,
                               residues = residues,
                               fraction = fraction)

ggplot(agg_polybasic_dt[datasource != "Kenyon"], aes(residues, fraction, fill = datasource)) + 
  geom_col(position = "dodge") + 
  scale_fill_manual(limits = c("Proteome", "Hartl"), 
                    labels = c("proteome", "agg"),
                    values = c("#999999", "#E31A1C"))
plot <- plot + theme_classic(20) + labs(y = "log2(frequency)", x = "Residue") +
  theme(legend.position = "none",
        axis.text = element_text(size = 16, color = "black"))

