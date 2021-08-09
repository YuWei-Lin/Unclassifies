### Compare the matlab result list to sil candidate GO terms
GOIDList <- read.csv("D:/GSE75688_Breast_Seed/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
DATA <- read.csv("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/GSE75688_BC_tumorCDF_ScorpVofZRM30.csv", header = T, stringsAsFactors = F)
DATA2 <- read.table("C:/Users/User/Desktop/MatLab Codes/breast_cancer_clusters_unique_GOterms.txt", sep = "\n", stringsAsFactors = F)
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

attach(DATA)
DATA <- DATA[order(DATA$Clust_p.Val, decreasing = T), ]
septerms <- DATA[ Clust_p.Val >= 0.95 , 1]
nsepterms <- DATA[ Clust_p.Val <= 0.05 , 1]

# Find the cluster number of where each septerm is located in Mat
Ovlapinsep <- 1:5 # Set up a dim number that can possibly accommondates the repeats of overlapping terms
for (a in 1:length(septerms)) {
  ECH <- which(sapply(Mat, function(x) septerms[a]%in%x))
  if(length(ECH) >= 1){
    for (b in 1:length(ECH)) {
      Tag <- septerms[a]
      if(as.numeric(ECH[b]) < 10){
        ECH[b] <- paste(0, ECH[b], sep = "")
      }
      Tag <- c(Tag, paste("cluster: ", ECH[b], sep = ""))
      Ovlapinsep <- data.frame.na(Ovlapinsep, Tag)
    }
  }
}
Ovlapinsep <- t(Ovlapinsep)
Ovlapinsep <- Ovlapinsep[2:nrow(Ovlapinsep), ]
Ovlapinsep <- Ovlapinsep[, !apply(is.na(Ovlapinsep), 2, all)] # remove all NA columns.
Ovlapinsep <- as.data.frame(Ovlapinsep, row.names = 1:nrow(Ovlapinsep))
Ovlapinsep <- Ovlapinsep[order(Ovlapinsep$V2), ]
NN <- GOIDList[GOIDList$ProcessName%in%Ovlapinsep$V1, ]
Ovlapinsep <- cbind(NN[match(Ovlapinsep$V1, NN$ProcessName), 2], Ovlapinsep)
Ovlapinsep <- cbind(Ovlapinsep, table(NN$GOID)[match(Ovlapinsep[,1], names(table(NN$GOID)))])
Ovlapinsep <- Ovlapinsep[,-1]
colnames(Ovlapinsep) <- c("GO terms", "K_Cluster", "GO-ID", "Gene Numbers")
write.table(Ovlapinsep, file = "GSE75688_Breast_Sil95ClusWalk_Comparison.csv", sep = ",", row.names = FALSE)

# Find the cluster number of where each nsepterm is located in Mat
Ovlapinnsep <- 1:5 # Set up a dim number that can possibly accommondates the repeats of overlapping terms
for (a in 1:length(nsepterms)) {
  ECH <- which(sapply(Mat, function(x) nsepterms[a]%in%x))
  if(length(ECH) >= 1){
    for (b in 1:length(ECH)) {
      Tag <- nsepterms[a]
      if(as.numeric(ECH[b]) < 10){
        ECH[b] <- paste(0, ECH[b], sep = "")
      }
      Tag <- c(Tag, paste("cluster: ", ECH[b], sep = ""))
      Ovlapinnsep <- data.frame.na(Ovlapinnsep, Tag)
    }
  }
}
Ovlapinnsep <- t(Ovlapinnsep)
Ovlapinnsep <- Ovlapinnsep[2:nrow(Ovlapinnsep), ]
Ovlapinnsep <- Ovlapinnsep[, !apply(is.na(Ovlapinnsep), 2, all)] # remove all NA columns.
Ovlapinnsep <- as.data.frame(Ovlapinnsep, row.names = 1:nrow(Ovlapinnsep))
Ovlapinnsep <- Ovlapinnsep[order(Ovlapinnsep$V2), ]
NN <- GOIDList[GOIDList$ProcessName%in%Ovlapinnsep$V1, ]
Ovlapinnsep <- cbind(NN[match(Ovlapinnsep$V1, NN$ProcessName), 2], Ovlapinnsep)
Ovlapinnsep <- cbind(Ovlapinnsep, table(NN$GOID)[match(Ovlapinnsep[,1], names(table(NN$GOID)))])
Ovlapinnsep <- Ovlapinnsep[,-1]
colnames(Ovlapinnsep) <- c("GO terms", "K_Cluster", "GO-ID", "Gene Numbers")
write.table(Ovlapinnsep, file = "GSE75688_Breast_Sil05ClusWalk_Comparison.csv", sep = ",", row.names = FALSE)

### Generate sil and tsne plots of the uniquely enriched clusters
Mel.stand <- read.csv("D:/GSE75688_Breast_Seed/GSE75688_Breast_NoBulk_CDF.csv", header = TRUE)
# Remove genes that contain zero more than 7/10
Mel.stand <- Mel.stand[apply(Mel.stand, 1, function(x) sum(is.na(x)) < (ncol(Mel.stand)*(0.7))), ]
Mel.stand <- t(Mel.stand)
# Molecular Subtypes labeling
library(Rtsne); library(gplots); library(cluster)
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe"
           , "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")
labels <- strsplit(rownames(Mel.stand), "_")
labels <- unlist(lapply(labels, function(x) x[1]))
tumors_colors=rep(colors[1], nrow(Mel.stand))
for(w in 2:length(table(labels))){
  tumors_colors[labels==names(table(labels))[w]]=colors[w]
}
# Cancer cells grouping
gps <- rep(1, nrow(Mel.stand))
for(g in 2:length(table(labels))){
  gps[labels==names(table(labels))[g]]=g
}
Mel.stand[is.na(Mel.stand)] <- 0

# Particularly generate tsne and sil plots for the clusters that contain uniquely enriched GO terms.
ClusN <- NULL
for (r in 1:78) {
  if(all(is.na(Mat[ , r])) == F){
    ClusN <- c(ClusN, paste("cluster ", r, sep = ""))
  }
}
ClusNN <- NULL
for (w in 1:length(ClusN)) {
  ID <- as.numeric(sub("cluster ", "", ClusN[w]))
  CHE <- sub(" ", "", strsplit(DATA2$V1[ID+1], ",")[[1]][2])
  X <- as.numeric(substr(CHE, as.numeric(gregexpr("K", CHE))+1, as.numeric(gregexpr("C", CHE))-1))
  Y <- as.numeric(substr(CHE, as.numeric(gregexpr("C", CHE))+1, nchar(CHE)))
  if( X < 10){
    X <- paste("0", X, sep = "")
  }
  if(Y < 10){
    Y <- paste("0", Y, sep = "")
  }
  CHE <- paste("K", X, "C", Y, sep = "")
  ClusNN <- c(ClusNN, CHE)
}

GNVa <- NULL
GNum <- NULL
Indsil <- NULL
GNAvgss <- NULL
for (r in 1:length(ClusNN)) {
  GENES <- read.csv(paste("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/ClusterGeneLists/GSE75688_Breast_ClusGenes_", ClusNN[r], ".csv", sep = ""), header = F, stringsAsFactors = F)
  GENES <- GENES$V1[3:nrow(GENES)]
  AD <- Mel.stand[ , colnames(Mel.stand)%in%GENES]
  AD <- AD[ ,!(apply(AD, 2, function(x) sum(!(is.na(x)))<=1|var(x ,na.rm = T)<0.0000000001))]
  AD <- AD[!(apply(AD, 1, function(x) sum(!(is.na(x)))<=1|var(x ,na.rm = T)<0.0000000001)), ]
  adjPP <- max((ncol(AD)%/%500)*5, 10)
  if(adjPP > 30){
    adjPP = 30
  }
  tsne <- Rtsne(AD, dims = 2, perplexity=adjPP, max_iter = 5000, check_duplicates = FALSE)
  cormat <- cor(t(AD), method = "pearson")
  cortrans <- 2*(1-cormat) 
  sil = silhouette (gps, cortrans)
  tt <- sil[ ,"sil_width"]
  GNsubsil <- NULL
  GNVa <- append(GNVa, length(tt[tt<=0])/length(tt), after = length(GNVa))
  for(k in 1:length(table(gps))){
    KK <- sil[ ,"cluster"] == k
    GNsubsil <- append(GNsubsil, sum(sil[KK,"sil_width"]<=0)/sum(KK), after = length(GNsubsil))
  }
  Indsil <- rbind(Indsil, GNsubsil)
  GNAvgss <- append(GNAvgss, mean(GNsubsil), after = length(GNAvgss))
  GNum <- append(GNum, ncol(AD), after = length(GNum))
    
  tiff(paste("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/Uniquely_Enriched_Clus_tsneplots/UCluster_", ClusNN[r],".tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
  plot(tsne$Y, main= paste("UCluster_", ClusNN[r], ncol(AD), format(length(tt[tt<=0])/length(tt), digits = 2, format = T)), cex.main = 0.8, col= tumors_colors, pch = 16, cex = 0.4)
  dev.off()
    
  tiff(paste("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/Uniquely_Enriched_Clus_Silplots/UCluster_", ClusNN[r],".tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
  plot(sil, main = paste("UCluster_", ClusNN[r], ncol(AD), format(length(tt[tt<=0])/length(tt), digits = 2, format = T)), col= colors[1:12], border=NA)
  dev.off()
}

# Randomly runs for computing p-value
p_time <- proc.time()
GNVaAcu <- NULL
GNVaAcuInd <- NULL
for(i in 1:length(GNum)){
  GNVaR <- NULL
  GNAvgssR <- NULL
  for(j in 1:1000){
    ADR <- Mel.stand[ ,sample(ncol(Mel.stand), GNum[i], replace = FALSE)]
    cormatR <- cor(t(ADR), method = "pearson")
    cortransR <- 2*(1-cormatR) 
    silR = silhouette (gps, cortransR)
    ttR <- silR[ ,"sil_width"]
    GNVaR <- append(GNVaR, length(ttR[ttR<=0])/length(ttR), after = length(GNVaR))
    GNsubsilR <- NULL
    for(k in 1:length(table(gps))){
      KKR <- silR[ ,"cluster"] == k
      GNsubsilR <- append(GNsubsilR, sum(silR[KKR,"sil_width"]<=0)/sum(KKR), after = length(GNsubsilR))
    }
    GNAvgssR <- append(GNAvgssR, mean(GNsubsilR), after = length(GNAvgssR))  
  }
  GNVaAcu <- append(GNVaAcu, 1-(sum(as.numeric(GNVa[i] > GNVaR))/1000), after = length(GNVaAcu))
  GNVaAcuInd <- append(GNVaAcuInd, 1-(sum(as.numeric(GNAvgss[i] > GNAvgssR))/1000), after = length(GNVaAcuInd))
}
t_time <- proc.time()-p_time
print(t_time)
E1 <- cbind(ClusN, ClusNN, GNum, Indsil)
E2 <- cbind(GNAvgss, GNVa, GNVaAcuInd, GNVaAcu)
E <- cbind(E1, E2)
E <- as.data.frame(E, stringsAsfactor = FALSE)
E <- E[rev(order(E$GNVaAcuInd)), ]
colnames(E) <- c("Cluster Index", "Enriched Cluster", "Gene Size", paste(names(table(labels))), "Clust_NSV", "Whole_NSV", "Clust_p-Val", "Whole_p-Val")
write.table(E, file = "GSE75688_BC_tumorCDF_ScorpVof_EnrichedCluster.csv", sep = ",", row.names = FALSE)


tiff(paste("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/Uniquely_Enriched_Clus_Silplots/UCluster_XX.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
plot(siltsne, main = paste("UCluster_XX", sep=""), col= colors[1:12], border=NA)
dev.off()