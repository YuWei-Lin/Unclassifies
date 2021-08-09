# Implement rejection sampling from an empirical distribution of 1D data.
# Inputs: 1D observation data, data value range, number of intervals in the value range, size of the random sample.
# Outputs: randomly sampled data drawn from the empirical distribution.
# Obtain the empirical PDF.
# Randomly draw numbers from a uniform distribution over the range.
# Accept a number a with probability Pr(X<=a) and Pr(X<=a) from empirical distribution.
# Continue until generating the sample of the designated size.

rejection_sampling <- function(empdatapoints, valrange, nintervals, nrsamplesize){
  # Calculate the PDF of empirical data.
  # Subdivide the valrange into nintervals intervals.  Calculate the PDF for each interval value.
  minx=min(valrange); maxx=max(valrange); dx=(maxx-minx)/nintervals
  intervals=seq(minx, maxx, dx)
  pdfvals=matrix(0, 1, nintervals+1)
  for (n in 1:nintervals) {
    val1=intervals[n]; val2=intervals[n+1]
    k=sum((empdatapoints>=val1)&(empdatapoints<=val2))
    pdfvals[n]=k/length(empdatapoints)
  }
  maxp = max(pdfvals)
  # Apply rejection sampling to generate random sampled points.
  rsampledpoints <- NULL
  while(length(rsampledpoints)<nrsamplesize){
    rvals=runif(nrsamplesize)
    rvals=(maxx-minx)*rvals+minx
    #rvals=rand(1,nrsamplesize)*dx+minx;
    qs=round((rvals-minx)/dx)+1
    qs[which(qs>(nintervals+1))]=nintervals+1
    p0s=pdfvals[qs] 
    p1s=runif(nrsamplesize)*maxp  
    ss=which(p1s<=p0s)
    rsampledpoints=c(rsampledpoints, rvals[ss])
  }
  rsampledpoints=rsampledpoints[1:nrsamplesize]
  # Handle the case when empdatapoints are concentrated in a value.
  if (length(unique(empdatapoints))==1){
    val=unique(empdatapoints)
    rsampledpoints=val*matrix(1, 1, nrsamplesize)
  }
  return(rsampledpoints)
}

# Evaluate the median score of a group of genes, and the significance of their distribution from a background distribution.
# Direction is determined by the mean score.
# Inputs: background scores, the data scores.
# Outputs: median score, significance of deviation.
# Define X as a random variable of correlation coefficients drawn from data.
# Define Y as a random variable of correlation coefficients drawn from the background distribution.
# Define Z=X-Y if X has positive deviation from background, and Z=Y-X vice versa.
# The significance score is Pr(Z>0)-Pr(Z<0).
# Difference from evaluate_score_significance.m: (1)Apply rejection sampling to draw random samples of the same size from X and Y distributions.  Evaluate Pr(Z>0) accordingly.

evaluate_score_significance2 <- function (bgscores, datascores, nintervals, nrsamplesize){
  n0=length(bgscores); n1=length(datascores)
  
  # Calculate the median score.
  X=datascores; medscore=median(X)
  
  # Calculate the significance of deviation.
  # Randomly sample bgscores ntrials times with the same dimension as datascores.
  # Skip when encounter invalid entries.
  
  if (medscore>=0){
    dir=1
  }
  else{
    dir=-1
  }
  
  # Generate nrsamplesize random samples from X and Y distributions.
  Y=bgscores
  minval=quantile(Y,0.05); maxval=quantile(Y,0.95)
  valrange=c(minval, maxval)
  rXs=rejection_sampling(X,valrange,nintervals,nrsamplesize)
  rYs=rejection_sampling(Y,valrange,nintervals,nrsamplesize)
  
  if (dir>0){
    val1=sum(rXs>rYs); val2=sum(rXs<rYs)
  }
  else if (dir<0){
    val1=sum(rXs<rYs); val2=sum(rXs>rYs)
  }
  else{
    val1=0; val2=0
  }
  pdiff=(val1-val2)/nrsamplesize
  return(c(medscore, pdiff))
}

### General Pipe needs a input case name
CaseN <- "GSE72056_Melanoma"
STAT <- c("Normal", "Cancerous")
REPL <- c("Mean", "Zero")
### File Directories
Filepath1 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", sep = "")
### Import data
for (Y in STAT) {
  if(Y == "Normal"){
    for (R in REPL){
      tar <- list.files(paste(Filepath1, Y, "/", sep = ""), pattern = "DEL")
      Mel.stand <- read.csv(paste(Filepath1, Y, "/", tar, sep = ""), header = T, stringsAsFactors = F)
      Mel.stand <- Mel.stand[!duplicated(Mel.stand[ , 1]), ]
      rownames(Mel.stand) <- Mel.stand$Gene
      Mel.stand <- Mel.stand[ ,-1] 
      Mel.stand <- Mel.stand[apply(Mel.stand, 1, function(x) sum(is.na(x)) < (ncol(Mel.stand)*(0.7))), ] #Remove genes that contain zero more than 7/10
      if(R == "Mean"){
        ### Convert NAs to "Mean"
        for (i in 1:nrow(Mel.stand)) { Mel.stand[ i, is.na(Mel.stand[i, ])] <- mean(na.omit(as.numeric(Mel.stand[i, ]))) }
        Filepath2 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", Y, "/K.Means/Filled with Mean/", sep = "")
        tar2 <- list.files(paste(Filepath2, "inputs/", sep = ""), pattern = "Membership")
        Mems <- read.csv(paste(Filepath2, "/inputs/", tar2, sep = ""), header = T, stringsAsFactors = F)
        tar3 <- list.files(paste(Filepath1, Y, "/", sep = ""), pattern = "AnnoRows")
        Anno <- read.csv(paste(Filepath1, Y, "/", tar3, sep = ""), header = T, stringsAsFactors = F)
        rownames(Anno) <- Anno$Sample
        Anno <- Anno[ ,-1]
        AnnoInd <- order(Anno[2, ], decreasing = F)
        Anno <- Anno[ , AnnoInd]
        Mel.stand <- Mel.stand[ ,AnnoInd]
        labels <- as.numeric(Anno[2, ]) 
        M_Index <- c(0, cumsum(table(labels)))
        Mel.stand <- t(Mel.stand)
        Medsc.Psigni <- NULL
        for (w in 2:ncol(Mems)) {
          AD <- Mel.stand[ , colnames(Mel.stand)%in%(Mems[ Mems[ ,w] == T, ]$Gene)]
          if(is.null(dim(AD)) == T){
            Medsc.Psigni <- rbind(Medsc.Psigni, c(NA, NA))
            next
          }
          #AD <- as.data.frame(AD)
          #rownames(AD) <- rownames(Mel.stand)
          #colnames(AD) <- Mems[ Mems[ ,w] == T, ]$Gene
          #AD <- AD[apply(AD, 1, function(x) sd(x)!=0), ]
          M_Index <- c(0, cumsum(table(labels)))
          cormat <- cor(t(AD), method = "pearson")
          cortrans <- 2*(1-cormat)
          cortrans[lower.tri(cortrans, diag = T)] <- NA
          intra.P.dist <- NULL
          inter.P.dist <- NULL
          for (p in 1:(length(M_Index)-1)) {
            SQ <- cortrans[(M_Index[p]+1):M_Index[p+1], (M_Index[p]+1):M_Index[p+1]]
            intra.P.dist <- c(intra.P.dist, as.vector(SQ[is.na(SQ)==F]))
            cortrans[(M_Index[p]+1):M_Index[p+1], (M_Index[p]+1):M_Index[p+1]] <- NA
          }
          inter.P.dist <- as.vector(cortrans[is.na(cortrans)==F])
          Medsc.Psigni <- rbind(Medsc.Psigni, evaluate_score_significance2(intra.P.dist, inter.P.dist, 100, 10000))
          print(w-1)
        }
        ClusIdx <- apply(Mems[ ,2:ncol(Mems)], 2, function(x) sum(x))
        Medsc.Psigni <- cbind(names(ClusIdx), ClusIdx, Medsc.Psigni)
        colnames(Medsc.Psigni) <- c("Cluster", "Gene Size", "interdis_Median", "pVal_Score")
        Medsc.Psigni <- Medsc.Psigni[order(Medsc.Psigni[ , 4], decreasing = T), ]
        write.csv(Medsc.Psigni, paste( Filepath2, CaseN, "_ClusHeterogeneityScore.csv", sep = ""), row.names = F)
        print(paste(CaseN, "___", Y, "_Cells_", R, "___Completed!"))
      }else{
        ### Convert NAs to "0"
        Mel.stand[is.na(Mel.stand)] = 0
        Filepath2 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", Y, "/K.Means/Filled with Zero/", sep = "")
        tar2 <- list.files(paste(Filepath2, "inputs/", sep = ""), pattern = "Membership")
        Mems <- read.csv(paste(Filepath2, "/inputs/", tar2, sep = ""), header = T, stringsAsFactors = F)
        tar3 <- list.files(paste(Filepath1, Y, "/", sep = ""), pattern = "AnnoRows")
        Anno <- read.csv(paste(Filepath1, Y, "/", tar3, sep = ""), header = T, stringsAsFactors = F)
        rownames(Anno) <- Anno$Sample
        Anno <- Anno[ ,-1]
        AnnoInd <- order(Anno[2, ], decreasing = F)
        Anno <- Anno[ , AnnoInd]
        Mel.stand <- Mel.stand[ ,AnnoInd]
        labels <- as.numeric(Anno[2, ]) 
        M_Index <- c(0, cumsum(table(labels)))
        Mel.stand <- t(Mel.stand)
        Medsc.Psigni <- NULL
        for (w in 2:ncol(Mems)) {
          AD <- Mel.stand[ , colnames(Mel.stand)%in%(Mems[ Mems[ ,w] == T, ]$Gene)]
          if(is.null(dim(AD)) == T){
            Medsc.Psigni <- rbind(Medsc.Psigni, c(NA, NA))
            next
          }
          #AD <- AD[apply(AD, 1, function(x) sd(x)!=0), ]
          M_Index <- c(0, cumsum(table(labels)))
          cormat <- cor(t(AD), method = "pearson")
          cortrans <- 2*(1-cormat)
          cortrans[lower.tri(cortrans, diag = T)] <- NA
          intra.P.dist <- NULL
          inter.P.dist <- NULL
          for (p in 1:(length(M_Index)-1)) {
            SQ <- cortrans[(M_Index[p]+1):M_Index[p+1], (M_Index[p]+1):M_Index[p+1]]
            intra.P.dist <- c(intra.P.dist, as.vector(SQ[is.na(SQ)==F]))
            cortrans[(M_Index[p]+1):M_Index[p+1], (M_Index[p]+1):M_Index[p+1]] <- NA
          }
          inter.P.dist <- as.vector(cortrans[is.na(cortrans)==F])
          Medsc.Psigni <- rbind(Medsc.Psigni, evaluate_score_significance2(intra.P.dist, inter.P.dist, 100, 10000))
          print(w-1)
        }
        ClusIdx <- apply(Mems[ ,2:ncol(Mems)], 2, function(x) sum(x))
        Medsc.Psigni <- cbind(names(ClusIdx), ClusIdx, Medsc.Psigni)
        colnames(Medsc.Psigni) <- c("Cluster", "Gene Size", "interdis_Median", "pVal_Score")
        Medsc.Psigni <- Medsc.Psigni[order(Medsc.Psigni[ , 4], decreasing = T), ]
        write.csv(Medsc.Psigni, paste( Filepath2, CaseN, "_ClusHeterogeneityScore.csv", sep = ""), row.names = F)
        print(paste(CaseN, "___", Y, "_Cells_", R, "___Completed!"))
      }
    }
  }else{
    for (R in REPL){
      tar <- list.files(paste(Filepath1, Y, "/", sep = ""), pattern = "DEL")
      Mel.stand <- read.csv(paste(Filepath1, Y, "/", tar, sep = ""), header = T, stringsAsFactors = F)
      Mel.stand <- Mel.stand[!duplicated(Mel.stand[ , 1]), ]
      rownames(Mel.stand) <- Mel.stand$Gene
      Mel.stand <- Mel.stand[ ,-1] 
      Mel.stand <- Mel.stand[apply(Mel.stand, 1, function(x) sum(is.na(x)) < (ncol(Mel.stand)*(0.7))), ] #Remove genes that contain zero more than 7/10
      if(R == "Mean"){
        ### Convert NAs to "Mean"
        for (i in 1:nrow(Mel.stand)) { Mel.stand[ i, is.na(Mel.stand[i, ])] <- mean(na.omit(as.numeric(Mel.stand[i, ]))) }
        Filepath2 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", Y, "/K.Means/Filled with Mean/", sep = "")
        tar2 <- list.files(paste(Filepath2, "inputs/", sep = ""), pattern = "Membership")
        Mems <- read.csv(paste(Filepath2, "/inputs/", tar2, sep = ""), header = T, stringsAsFactors = F)
        tar3 <- list.files(paste(Filepath1, Y, "/", sep = ""), pattern = "AnnoRows")
        Anno <- read.csv(paste(Filepath1, Y, "/", tar3, sep = ""), header = T, stringsAsFactors = F)
        rownames(Anno) <- Anno$Sample
        Anno <- Anno[ ,-1]
        AnnoInd <- order(Anno[1, ], decreasing = F)
        Anno <- Anno[ , AnnoInd]
        Mel.stand <- Mel.stand[ ,AnnoInd]
        labels <- as.numeric(Anno[1, ]) 
        M_Index <- c(0, cumsum(table(labels)))
        Mel.stand <- t(Mel.stand)
        Medsc.Psigni <- NULL
        for (w in 2:ncol(Mems)) {
          AD <- Mel.stand[ , colnames(Mel.stand)%in%(Mems[ Mems[ ,w] == T, ]$Gene)]
          if(is.null(dim(AD)) == T){
            Medsc.Psigni <- rbind(Medsc.Psigni, c(NA, NA))
            next
          }
          #AD <- AD[apply(AD, 1, function(x) sd(x)!=0), ]
          M_Index <- c(0, cumsum(table(labels)))
          cormat <- cor(t(AD), method = "pearson")
          cortrans <- 2*(1-cormat)
          cortrans[lower.tri(cortrans, diag = T)] <- NA
          intra.P.dist <- NULL
          inter.P.dist <- NULL
          for (p in 1:(length(M_Index)-1)) {
            SQ <- cortrans[(M_Index[p]+1):M_Index[p+1], (M_Index[p]+1):M_Index[p+1]]
            intra.P.dist <- c(intra.P.dist, as.vector(SQ[is.na(SQ)==F]))
            cortrans[(M_Index[p]+1):M_Index[p+1], (M_Index[p]+1):M_Index[p+1]] <- NA
          }
          inter.P.dist <- as.vector(cortrans[is.na(cortrans)==F])
          Medsc.Psigni <- rbind(Medsc.Psigni, evaluate_score_significance2(intra.P.dist, inter.P.dist, 100, 10000))
          print(w-1)
        }
        ClusIdx <- apply(Mems[ ,2:ncol(Mems)], 2, function(x) sum(x))
        Medsc.Psigni <- cbind(names(ClusIdx), ClusIdx, Medsc.Psigni)
        colnames(Medsc.Psigni) <- c("Cluster", "Gene Size", "interdis_Median", "pVal_Score")
        Medsc.Psigni <- Medsc.Psigni[order(Medsc.Psigni[ , 4], decreasing = T), ]
        write.csv(Medsc.Psigni, paste( Filepath2, CaseN, "_ClusHeterogeneityScore.csv", sep = ""), row.names = F)
        print(paste(CaseN, "___", Y, "_Cells_", R, "___Completed!"))
      }else{
        ### Convert NAs to "0"
        Mel.stand[is.na(Mel.stand)] = 0
        Filepath2 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", Y, "/K.Means/Filled with Zero/", sep = "")
        tar2 <- list.files(paste(Filepath2, "inputs/", sep = ""), pattern = "Membership")
        Mems <- read.csv(paste(Filepath2, "/inputs/", tar2, sep = ""), header = T, stringsAsFactors = F)
        tar3 <- list.files(paste(Filepath1, Y, "/", sep = ""), pattern = "AnnoRows")
        Anno <- read.csv(paste(Filepath1, Y, "/", tar3, sep = ""), header = T, stringsAsFactors = F)
        rownames(Anno) <- Anno$Sample
        Anno <- Anno[ ,-1]
        AnnoInd <- order(Anno[1, ], decreasing = F)
        Anno <- Anno[ , AnnoInd]
        Mel.stand <- Mel.stand[ ,AnnoInd]
        labels <- as.numeric(Anno[1, ]) 
        M_Index <- c(0, cumsum(table(labels)))
        Mel.stand <- t(Mel.stand)
        Medsc.Psigni <- NULL
        for (w in 2:ncol(Mems)) {
          AD <- Mel.stand[ , colnames(Mel.stand)%in%(Mems[ Mems[ ,w] == T, ]$Gene)]
          if(is.null(dim(AD)) == T){
            Medsc.Psigni <- rbind(Medsc.Psigni, c(NA, NA))
            next
          }
          #AD <- AD[apply(AD, 1, function(x) sd(x)!=0), ]
          M_Index <- c(0, cumsum(table(labels)))
          cormat <- cor(t(AD), method = "pearson")
          cortrans <- 2*(1-cormat)
          cortrans[lower.tri(cortrans, diag = T)] <- NA
          intra.P.dist <- NULL
          inter.P.dist <- NULL
          for (p in 1:(length(M_Index)-1)) {
            SQ <- cortrans[(M_Index[p]+1):M_Index[p+1], (M_Index[p]+1):M_Index[p+1]]
            intra.P.dist <- c(intra.P.dist, as.vector(SQ[is.na(SQ)==F]))
            cortrans[(M_Index[p]+1):M_Index[p+1], (M_Index[p]+1):M_Index[p+1]] <- NA
          }
          inter.P.dist <- as.vector(cortrans[is.na(cortrans)==F])
          Medsc.Psigni <- rbind(Medsc.Psigni, evaluate_score_significance2(intra.P.dist, inter.P.dist, 100, 10000))
          print(w-1)
        }
        ClusIdx <- apply(Mems[ ,2:ncol(Mems)], 2, function(x) sum(x))
        Medsc.Psigni <- cbind(names(ClusIdx), ClusIdx, Medsc.Psigni)
        colnames(Medsc.Psigni) <- c("Cluster", "Gene Size", "interdis_Median", "pVal_Score")
        Medsc.Psigni <- Medsc.Psigni[order(Medsc.Psigni[ , 4], decreasing = T), ]
        write.csv(Medsc.Psigni, paste( Filepath2, CaseN, "_ClusHeterogeneityScore.csv", sep = ""), row.names = F)
        print(paste(CaseN, "___", Y, "_Cells_", R, "___Completed!"))
      }
    }
  }
}


