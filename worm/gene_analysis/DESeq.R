BiocManager::install("DESeq2")
library(data.table)
library(ggplot2)
library(DESeq2)
library(plyr)
library(dplyr)


### Calculate TPM ###
N2_dtA <- readRDS("N2_dtA.rds")
temp <- N2_dtA[position > 15 & stopdist < -15, c(1:14)]
N2samplesA <- names(N2_dtA[, c(8:14)])
for(i in N2samplesA) {
  print(i)
  temp[, paste0(i,"_sum") := lapply(.SD[, ..i], sum), by = orf]
  temp[, length1 := (length - 30) * 3 / 1000]
  temp[, rpk1 := temp[[paste0(i,"_sum")]] / length1]
  temp[, paste0(i,"_tpm") := (rpk1 / sum(temp[position == 16]$rpk1)) * 10^6]
}
N2_expression <- temp[, .SD[which.min(position)], by = orf]
#saveRDS(N2_expression, "N2_expression.rds")


### Test correlation ###
ggplot(N2_expression, aes(D1_1A_tpm, D1_2A_tpm)) + geom_point() +
  scale_x_log10(limits = c(1e-5,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-5,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
  annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))

cor(N2_expression[D1_1A_sum > 64 & D1_2A_sum > 64]$D1_1A_tpm,
    N2_expression[D1_1A_sum > 64 & D1_2A_sum > 64]$D1_2A_tpm, method = "pearson")
cor(N2_expression[D12_1A_sum > 64 & D12_2A_sum > 64]$D12_1A_tpm,
    N2_expression[D12_1A_sum > 64 & D12_2A_sum > 64]$D12_2A_tpm, method = "pearson")

cor(N2_expression[D1_1A_sum > 64 & D12_1A_sum > 64]$D1_1A_tpm,
    N2_expression[D1_1A_sum > 64 & D12_1A_sum > 64]$D12_1A_tpm, method = "pearson")
cor(N2_expression[D1_2A_sum > 64 & D12_2A_sum > 64]$D1_2A_tpm,
    N2_expression[D1_2A_sum > 64 & D12_2A_sum > 64]$D12_2A_tpm, method = "pearson")


### DESeq ###
# Create a matrix with the read counts
N2_expression <- readRDS("N2_expression.rds")
rawcountData <- as.matrix(N2_expression[, .(D12_1A_sum, D12_2A_sum, D1_1A_sum, D1_2A_sum)])
row.names(rawcountData) <- as.character(N2_expression$orf)
rawcountData = rawcountData[rowSums(rawcountData) > 1, ]

# Define your variables
condition <- factor(c('D12','D12','D1','D1'))

# Create the colData object
colData <- data.frame(condition = condition)
row.names(colData) <- c('D12_1A_sum', 'D12_2A_sum', 'D1_1A_sum', 'D1_2A_sum')
  
# Create DEseqset
dds <- DESeqDataSetFromMatrix(countData = rawcountData,
                              colData = colData,
                              design = ~ condition)
  
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c("condition", "D12", "D1"))


# This block will return a data.table with DESeq output
deseq_gene <- as.data.frame(res)
setDT(deseq_gene, keep.rownames = T)
colnames(deseq_gene)[1] <- "orf"
setkeyv(deseq_gene, c("orf"))
summary(res)
View(deseq_gene[log2FoldChange > 1 & log2FoldChange < Inf & padj < 0.05])

N2_expression_deseq <- N2_expression[deseq_gene]
N2_expression_deseq[, group := ifelse(padj < 0.05 & log2FoldChange >= 1, 1, 0)]
N2_expression_deseq[, group := ifelse(padj < 0.05 & log2FoldChange <= -1, 2, group)]
N2_expression_deseq[, D1_tpm := (D1_1A_tpm + D1_2A_tpm) / 2]
N2_expression_deseq[, D12_tpm := (D12_1A_tpm + D12_2A_tpm) / 2]
genenames <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/geneNames.csv', header = TRUE, stringsAsFactors = TRUE))
i <- cbind(match(N2_expression_deseq$orf, genenames$Transcript.stable.ID))
N2_expression_deseq <- cbind(N2_expression_deseq, WBgene = genenames[i]$Gene.stable.ID)
hartl_proteome <- as.data.table(read.csv('/Users/KevinStein/Desktop/Lab/Lab_notebook/1_Aging/Worm/Expt01_MassSpec analysis of age-dependent aggregates/Analysis/proteinLists/hartl_proteome.csv', header = TRUE, stringsAsFactors = TRUE))
i <- cbind(match(N2_expression_deseq$WBgene, hartl_proteome$WBgene))
N2_expression_deseq <- cbind(N2_expression_deseq, gene = hartl_proteome[i]$gene)
N2_expression_deseq <- cbind(N2_expression_deseq, uniprot = hartl_proteome[i]$uniprot)
#saveRDS(N2_expression_deseq, "N2_expression_deseq.rds")

write.csv(N2_expression_deseq[D1_1A_sum >= 64 & D1_2A_sum >= 64 & 
                                D12_1A_sum >= 64 & D12_2A_sum >= 64 & 
                                log2FoldChange < Inf & log2FoldChange > -Inf & padj < 0.05], "N2_expression.csv")

N2_expression_deseq <- readRDS("N2_expression_deseq.rds")
