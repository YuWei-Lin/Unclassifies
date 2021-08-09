
GSERAW <- read.table(file.choose(), fill = TRUE, header = TRUE, quote = "\"", sep = "\t", stringsAsFactors = FALSE)

GSEArray <- read.table(file.choose(), fill = TRUE, header = TRUE, quote = "\"", sep = "\t")

sort(table(unlist(lapply(strsplit(colnames(GSERAW), "\\..."), function(x) x[1]))))

#strsplit(gsub("\\."," ","Signomod.colon...5643"),"   ")

Bre <- GSERAW[ ,grepl("Breast", colnames(GSERAW))]
Bre <- GSEArray[ ,colnames(GSEArray)%in%Bre[1,]]
Bre <- cbind(GSEArray$ID_REF, Bre)
Col <- GSERAW[ ,grepl("Colon", colnames(GSERAW))|grepl("colon", colnames(GSERAW))]
Col <- GSEArray[ ,colnames(GSEArray)%in%Col[1,]]
Col <- cbind(GSEArray$ID_REF, Col)
Kid <- GSERAW[ ,grepl("Kidney", colnames(GSERAW))]
Kid <- GSEArray[ ,colnames(GSEArray)%in%Kid[1,]]
Kid <- cbind(GSEArray$ID_REF, Kid)
Ova <- GSERAW[ ,grepl("Ovary", colnames(GSERAW))]
Ova <- GSEArray[ ,colnames(GSEArray)%in%Ova[1,]]
Ova <- cbind(GSEArray$ID_REF, Ova)
Ute <- GSERAW[ ,grepl("Uterus", colnames(GSERAW))]
Ute <- GSEArray[ ,colnames(GSEArray)%in%Ute[1,]]
Ute <- cbind(GSEArray$ID_REF, Ute)
Lun <- GSERAW[ ,grepl("Lung", colnames(GSERAW))]
Lun <- GSEArray[ ,colnames(GSEArray)%in%Lun[1,]]
Lun <- cbind(GSEArray$ID_REF, Lun)
Pro <- GSERAW[ ,grepl("Prostate", colnames(GSERAW))]
Pro <- GSEArray[ ,colnames(GSEArray)%in%Pro[1,]]
Pro <- cbind(GSEArray$ID_REF, Pro)
Ome <- GSERAW[ ,grepl("Omentum", colnames(GSERAW))]
Ome <- GSEArray[ ,colnames(GSEArray)%in%Ome[1,]]
Ome <- cbind(GSEArray$ID_REF, Ome)
End <- GSERAW[ ,grepl("Endometrium", colnames(GSERAW))]
End <- GSEArray[ ,colnames(GSEArray)%in%End[1,]]
End <- cbind(GSEArray$ID_REF, End)
Liv <- GSERAW[ ,grepl("Liver", colnames(GSERAW))]
Liv <- GSEArray[ ,colnames(GSEArray)%in%Liv[1,]]
Liv <- cbind(GSEArray$ID_REF, Liv)
Rec <- GSERAW[ ,grepl("Rectum", colnames(GSERAW))]
Rec <- GSEArray[ ,colnames(GSEArray)%in%Rec[1,]]
Rec <- cbind(GSEArray$ID_REF, Rec)
Thy <- GSERAW[ ,grepl("Thyroid", colnames(GSERAW))]
Thy <- GSEArray[ ,colnames(GSEArray)%in%Thy[1,]]
Thy <- cbind(GSEArray$ID_REF, Thy)
Cer <- GSERAW[ ,grepl("Cervix", colnames(GSERAW))]
Cer <- GSEArray[ ,colnames(GSEArray)%in%Cer[1,]]
Cer <- cbind(GSEArray$ID_REF, Cer)

write.table(Cer, file = "GSE2109Info2-Cervix.txt", row.names = FALSE, sep = "\t")
write.table(Cer, file = "GSE2109Final-Cervix.txt", row.names = FALSE, sep = "\t")


