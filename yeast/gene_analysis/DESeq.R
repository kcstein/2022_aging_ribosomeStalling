source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(DESeq2)


### Calculate TPM ###
sc_dtA <- readRDS("sc_dtA.rds")
temp <- sc_dtA[position > 15 & stopdist < -15, c(1:19)]
samplesA <- names(sc_dtA[, c(8:19)])
for(i in samplesA) {
  print(i)
  temp[, paste0(i,"_sum") := lapply(.SD[, ..i], sum), by = orf]
  temp[, length1 := (length - 30) * 3 / 1000]
  temp[, rpk1 := temp[[paste0(i,"_sum")]] / length1]
  temp[, paste0(i,"_tpm") := (rpk1 / sum(temp[position == 16]$rpk1)) * 10^6]
}
expression <- temp[, .SD[which.min(position)], by = orf]
#saveRDS(expression, "expression.rds")


### Test correlation ###
ggplot(expression, aes(WT_D0_tpm, WT_D4_tpm)) + geom_point() +
  scale_x_log10(limits = c(1e-3,1e5), breaks = c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-3","-2","-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-3,1e5), breaks = c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-3","-2","-1", "0","1", "2", "3","4","5")) + 
  annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))

# Comparing replicates
cor(expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64]$WT_D0_1A_tpm,
    expression[WT_D0_1A_sum > 64 & WT_D0_2A_sum > 64]$WT_D0_2A_tpm, method = "pearson")
cor(expression[WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64]$WT_D4_1A_tpm,
    expression[WT_D4_1A_sum > 64 & WT_D4_2A_sum > 64]$WT_D4_2A_tpm, method = "pearson")
cor(expression[sch9_D0_1A_sum > 64 & sch9_D0_2A_sum > 64]$sch9_D0_1A_tpm,
    expression[sch9_D0_1A_sum > 64 & sch9_D0_2A_sum > 64]$sch9_D0_2A_tpm, method = "pearson")
cor(expression[sch9_D4_1A_sum > 64 & sch9_D4_2A_sum > 64]$sch9_D4_1A_tpm,
    expression[sch9_D4_1A_sum > 64 & sch9_D4_2A_sum > 64]$sch9_D4_2A_tpm, method = "pearson")

# Comparing aged to young
cor(expression[WT_D0_sum > 128 & WT_D4_sum > 128]$WT_D0_tpm,
    expression[WT_D0_sum > 128 & WT_D4_sum > 128]$WT_D4_tpm, method = "pearson")
cor(expression[sch9_D0_sum > 128 & sch9_D4_sum > 128]$sch9_D0_tpm,
    expression[sch9_D0_sum > 128 & sch9_D4_sum > 128]$sch9_D4_tpm, method = "pearson")

# Comparing WT to sch9
cor(expression[WT_D0_sum > 128 & sch9_D0_sum > 128]$WT_D0_tpm,
    expression[WT_D0_sum > 128 & sch9_D0_sum > 128]$sch9_D0_tpm, method = "pearson")
cor(expression[WT_D4_sum > 128 & sch9_D4_sum > 128]$WT_D4_tpm,
    expression[WT_D4_sum > 128 & sch9_D4_sum > 128]$sch9_D4_tpm, method = "pearson")


### DESeq ###
# Create a matrix with the read counts
expression <- readRDS("expression.rds")
rawcountData <- as.matrix(expression[, .(WT_D4_1A_sum, WT_D4_2A_sum, WT_D0_1A_sum, WT_D0_2A_sum)])
row.names(rawcountData) <- as.character(expression$orf)
rawcountData = rawcountData[rowSums(rawcountData) > 1, ]

# Define your variables
condition <- factor(c('D4','D4','D0','D0'))

# Create the colData object
colData <- data.frame(condition = condition)
row.names(colData) <- c('WT_D4_1A_sum', 'WT_D4_2A_sum', 'WT_D0_1A_sum', 'WT_D0_2A_sum')
  
# Create DEseqset
dds <- DESeqDataSetFromMatrix(countData = rawcountData,
                              colData = colData,
                              design = ~ condition)
  
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c("condition", "D4", "D0"))


# This block will return a data.table with DESeq output
deseq_gene <- as.data.frame(res)
setDT(deseq_gene, keep.rownames = T)
colnames(deseq_gene)[1] <- "orf"
setkeyv(deseq_gene, c("orf"))
summary(res)
View(deseq_gene[log2FoldChange > 1 & log2FoldChange < Inf & padj < 0.05])
WT_expression <- expression[deseq_gene]
WT_expression[, group := ifelse(padj < 0.05 & log2FoldChange >= 1, 1, 0)]
WT_expression[, group := ifelse(padj < 0.05 & log2FoldChange <= -1, 2, group)]
WT_expression[, D0_tpm := (WT_D0_1A_tpm + WT_D0_2A_tpm) / 2]
WT_expression[, D4_tpm := (WT_D4_1A_tpm + WT_D4_2A_tpm) / 2]
#saveRDS(WT_expression, "WT_expression.rds")

write.csv(WT_expression[WT_D0_1A_sum >= 64 & WT_D0_2A_sum >= 64 & 
                       WT_D4_1A_sum >= 64 & WT_D4_2A_sum >= 64 & 
                       log2FoldChange < Inf & log2FoldChange > -Inf & padj < 0.05], "expression.csv")

WT_expression <- readRDS("WT_expression.rds")




### Sch9 DESeq ###
# Create a matrix with the read counts
expression <- readRDS("expression.rds")
rawcountData_sch9 <- as.matrix(expression[, .(sch9_D4_1A_sum, sch9_D4_2A_sum, sch9_D0_1A_sum, sch9_D0_2A_sum)])
row.names(rawcountData_sch9) <- as.character(expression$orf)
rawcountData_sch9 = rawcountData_sch9[rowSums(rawcountData_sch9) > 1, ]
condition <- factor(c('D4','D4','D0','D0'))
colData_sch9 <- data.frame(condition = condition)
row.names(colData_sch9) <- c('sch9_D4_1A_sum', 'sch9_D4_2A_sum', 'sch9_D0_1A_sum', 'sch9_D0_2A_sum')

# Create DEseqset
dds_sch9 <- DESeqDataSetFromMatrix(countData = rawcountData_sch9,
                              colData = colData_sch9,
                              design = ~ condition)
dds_sch9 <- DESeq(dds_sch9)
res_sch9 <- results(dds_sch9, alpha = 0.05, contrast = c("condition", "D4", "D0"))


# This block will return a data.table with DESeq output
deseq_gene_sch9 <- as.data.frame(res_sch9)
setDT(deseq_gene_sch9, keep.rownames = T)
colnames(deseq_gene_sch9)[1] <- "orf"
setkeyv(deseq_gene_sch9, c("orf"))
summary(res_sch9)
View(deseq_gene_sch9[log2FoldChange > 1 & log2FoldChange < Inf & padj < 0.05])
sch9_expression <- expression[deseq_gene_sch9]
saveRDS(sch9_expression, "sch9_expression.rds")
# write.csv(sch9_expression[sch9_D0_1A_sum >= 64 & sch9_D0_2A_sum >= 64 & 
#                             sch9_D4_1A_sum >= 64 & sch9_D4_2A_sum >= 64 & 
#                        log2FoldChange < Inf & log2FoldChange > -Inf & padj < 0.05], "sch9_expression.csv")
write.csv(sch9_expression, "sch9_expressionAll.csv")

