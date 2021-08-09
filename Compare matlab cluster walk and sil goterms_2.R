### Compare the matlab result list to sil candidate GO terms
GOIDList <- read.csv("D:/GSE75688_Breast_Seed/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
DATA <- read.csv("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/Inputs/GSE75688_Breast_KClusters_GOWhole.csv", header = T, stringsAsFactors = F)
DATA2 <- read.table("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/breast_cancer_clusters_unique_GOterms.txt", sep = "\n", stringsAsFactors = F)
DATA3 <- DATA2[ -1:-986, ] # Need to identify the row number of first cluster involving enriched GOterms in matlab output txt.
PosOfenriClus <- which(lapply(strsplit(DATA3, ' '), function(x) which(x == 'cluster'))==1)
PosOfenriClus <- c(PosOfenriClus, length(DATA3)+1)
Mat <- 1:max(diff(PosOfenriClus))
for (i in 1:length(PosOfenriClus)) {
  Mat2 <- DATA3[PosOfenriClus[i]:(PosOfenriClus[i+1]-1)]
  Mat <- data.frame.na(Mat, Mat2)
} # Here error is negligible!
Mat <- Mat[ ,-1]
colnames(Mat) <- Mat[1, ]
Mat <- Mat[-1, ]
rownames(Mat) <- 1:nrow(Mat)
for (i in 1:nrow(Mat)) {
  for (j in 1:ncol(Mat)) {
    if(is.na(Mat[i,j]) == F){
      Mat[i,j] <- gsub(")", "", strsplit(Mat[i,j], " ")[[1]][4])
    }
  }
}
GOIDList <- unique.data.frame(GOIDList[ ,-1])

### Change GO categories to GOID
UniClustVal <- NULL
for (w in 1:78) {
  print(w)
  Uni.FDR <- DATA[ DATA$Enrichment.GOTerms%in%Mat[ , w], w+2]
  Block.Clus <- cbind(GOIDList[GOIDList$ProcessName%in%Mat[ ,w], 1], Uni.FDR)
  UniClustVal <- cbind.na(UniClustVal, Block.Clus)
}

UniClustVal <- UniClustVal[ ,-1]
colnames(UniClustVal)[seq(1:ncol(UniClustVal))%%2!=0] <- colnames(Mat)
write.csv(UniClustVal, "C:/Users/User/Desktop/20190402 Progress/Breast_GO.REVIGO.csv", row.names = FALSE)
