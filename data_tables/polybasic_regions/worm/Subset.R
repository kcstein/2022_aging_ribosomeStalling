library(data.table)


dir <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/"
files <- list.files(dir, pattern = "*.csv")
samples <- sub("\\.csv", '', files)

for(i in 1:length(files)) {
  print(files[i])
  dt <- as.data.table(read.csv(paste0(dir, files[i]), header = FALSE, stringsAsFactors = TRUE))
  colnames(dt) <- c("orf", "position")
  dt[, name := paste0(orf, "_", position)]
  assign(samples[i], dt)
}

KR5of10 <- KR5of10[!KR5of10$name %in% KR6of10$name]
KR4of10 <- KR4of10[!KR4of10$name %in% KR6of10$name]
KR4of10 <- KR4of10[!KR4of10$name %in% KR5of10$name]
KR3of10 <- KR3of10[!KR3of10$name %in% KR6of10$name]
KR3of10 <- KR3of10[!KR3of10$name %in% KR5of10$name]
KR3of10 <- KR3of10[!KR3of10$name %in% KR4of10$name]

KR5of5 <- KR5of5[!KR5of5$name %in% KR6of6$name]
KR4of4 <- KR4of4[!KR4of4$name %in% KR5of5$name]
KR4of4 <- KR4of4[!KR4of4$name %in% KR6of6$name]
KR3of3 <- KR3of3[!KR3of3$name %in% KR6of6$name]
KR3of3 <- KR3of3[!KR3of3$name %in% KR5of5$name]
KR3of3 <- KR3of3[!KR3of3$name %in% KR4of4$name]

dfs <- Filter(function(x) is(x, "data.table"), mget(ls()))

for(i in 1:length(dfs)) {
  write.csv(dfs[i], paste0(dir, names(dfs[i]), "_only.csv"), row.names = FALSE)
}

DE4of4 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/DE4of4.csv", stringsAsFactors = TRUE, header = FALSE))
DE5of5 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/DE5of5.csv", stringsAsFactors = TRUE, header = FALSE))
DE5of10 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/DE5of10.csv", stringsAsFactors = TRUE, header = FALSE))
DE6of10 <- as.data.table(read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Worm/PolybasicRegions/DE6of10.csv", stringsAsFactors = TRUE, header = FALSE))
setnames(DE4of4, c("orf", "position"))
setnames(DE5of5, c("orf", "position"))
setnames(DE5of10, c("orf", "position"))
setnames(DE6of10, c("orf", "position"))
DE4of4[, name := paste0(orf, "_", position)]
DE5of5[, name := paste0(orf, "_", position)]
DE5of10[, name := paste0(orf, "_", position)]
DE6of10[, name := paste0(orf, "_", position)]
DE4of4 <- DE4of4[!DE4of4$name %in% DE5of5$name]
DE5of10 <- DE5of10[!DE5of10$name %in% DE6of10$name]
write.csv(DE4of4, paste0(dir, "DE4of4_only.csv"), row.names = FALSE)
write.csv(DE5of5, paste0(dir, "DE5of5_only.csv"), row.names = FALSE)
write.csv(DE5of10, paste0(dir, "DE5of10_only.csv"), row.names = FALSE)
write.csv(DE6of10, paste0(dir, "DE6of10_only.csv"), row.names = FALSE)
