library(data.table)


dir <- "/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Yeast/PolybasicRegions/"
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
KR4of4 <- KR4of4[!KR4of4$name %in% KR6of6$name]
KR4of4 <- KR4of4[!KR4of4$name %in% KR5of5$name]
KR3of3 <- KR3of3[!KR3of3$name %in% KR6of6$name]
KR3of3 <- KR3of3[!KR3of3$name %in% KR5of5$name]
KR3of3 <- KR3of3[!KR3of3$name %in% KR4of4$name]

rm(dt)
dfs <- Filter(function(x) is(x, "data.table"), mget(ls()))

for(i in 1:length(dfs)) {
  write.csv(dfs[i], paste0(dir, names(dfs[i]), "_only.csv"), row.names = FALSE)
}

