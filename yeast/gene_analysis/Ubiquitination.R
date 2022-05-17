### Correlation between stalling and ubiquitination ###
WT_stalling_peaks_ageDep <- readRDS("doc/WT_stalling_peaks_ageDep.rds")
ubiq <- as.data.table(read.csv("/Users/KevinStein/Desktop/Manuscript/Data/datasets/Duttler.csv", header = TRUE, stringsAsFactors = TRUE))
setkeyv(ubiq, "YORF")
proteome_orfs <- as.data.table(unique(sc_fishers$orf))
setnames(proteome_orfs, "orf")
ubiq_subset <- ubiq[ubiq$YORF %in% proteome_orfs$orf] # restricting ubiquitin dataset to the orfs that had stalling Fisher test applied
stall <- as.data.table(unique(WT_stalling_peaks_ageDep$orf))
setnames(stall, "orf")

length(ubiq_subset[ubiq_subset$YORF %in% stall$orf]$YORF) # 140 (stall and ubiq)
length(stall[!stall$orf %in% ubiq_subset$YORF]$orf) # 797 (stall, no ubiq)
length(ubiq_subset[!ubiq_subset$YORF %in% stall$orf]$YORF) # 351 (ubiq, no stall)
noUbiq <- proteome_orfs[!proteome_orfs$orf %in% ubiq_subset$YORF]
length(noUbiq[!noUbiq$orf %in% stall$orf]$orf) # 2794 (no stall, no ubiq)

counts <- matrix(c(140, 797, 351, 2794), nrow = 2)
fisher.test(counts)

# fraction of proteome that is ubiquitinated but not stall
# fraction of stalling that is ubiquitinated
temp <- data.table(Category = c("No", "Stall"), 
                   FractionUbiq = c(length(ubiq_subset[!ubiq_subset$YORF %in% stall$orf]$YORF) /
                                      length(proteome_orfs[!proteome_orfs$orf %in% stall$orf]$orf),
                                    length(ubiq_subset[ubiq_subset$YORF %in% stall$orf]$YORF) /
                                      length(stall$orf)))

plot <- ggplot(temp, aes(x = Category, y=FractionUbiq, fill = Category)) + geom_col(color = "black", size = 0.5) + scale_y_continuous(limits = c(0,0.15), expand = c(0.0008,0.0008)) +
  scale_fill_manual(limits = c("No","Stall"), values = c("#80B1D3", "#FDB462"), name = "")
plot <- plot + theme_classic(18) + labs(y = "Fraction of ubiquitinated\nproteins in dataset", x = "") +
  theme(legend.position = "none", legend.background = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Stalling_DuttlerUbiq_0.002388.pdf", plot, width = 3.2, height = 4, dpi = 300, useDingbats = F)
