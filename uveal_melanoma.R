#untar('gdac.broadinstitute.org_UVM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_Clinical.Level_1.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM-TP.Aggregate_AnalysisFeatures.Level_4.2015082100.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz', exdir = '.')

data <- read.delim('./gdac.broadinstitute.org_UVM-TP.Aggregate_AnalysisFeatures.Level_4.2015082100.0.0/UVM-TP.transferedsamplefeatures.txt')

data[1:5,1:5]
data_gdac <- t(data[,-1])
colnames(data_gdac) <- data[,1]
rownames(data_gdac) <- colnames(data)[-1]
data <- data.frame(data_gdac[-dim(data_gdac)[1],])

data[, c(1:4, 15:250)] <- sapply(data[, c(1:4, 15:250)], as.character)
data[, c(1:4, 15:250)] <- sapply(data[, c(1:4, 15:250)], as.numeric)
data$ID <- rownames(data)

summary(data)
colnames(data)

ID_class1 <- data$ID[data$CLUS_CN_cNMF == 1]
ID_class2 <- data$ID[data$CLUS_CN_cNMF == 2]
ID_class3 <- data$ID[data$CLUS_CN_cNMF == 3]

boxplot(data$CN_3p_Del ~ as.factor(data$CLUS_CN_cNMF))


########### Cluster CN NMF ##########################################

################# RNA for GSEA ######################################
############http://www.hammerlab.org/################################
library(RTCGAToolbox)

getFirehoseDatasets()
getFirehoseRunningDates(last = NULL)
readData = getFirehoseData (dataset="UVM", runDate="20151101",forceDownload = TRUE,
                            Clinic=TRUE, Mutation=FALSE, Methylation=FALSE, RNAseq2_Gene_Norm=TRUE)

clin = getData(readData, "Clinical")
mRNA <- getData(readData, 'RNASeq2GeneNorm')

colnames(mRNA) <- gsub('-...-...-....-..', '', colnames(mRNA))
colnames(mRNA) <- gsub('-', '\\.', colnames(mRNA))
colnames(mRNA)

data_GSEA <- c()
data_GSEA <- cbind(mRNA[, colnames(mRNA) %in% ID_class1],
                   mRNA[, colnames(mRNA) %in% ID_class2],
                   mRNA[, colnames(mRNA) %in% ID_class3])

GSEA <- data.frame(NAME = rownames(data_GSEA), 
                   DESCRIPTION = rep('na', dim(data_GSEA)[1]), data_GSEA)
cls <- c(rep(0, length(ID_class1)),
         rep(1, length(ID_class2)),
         rep(2, length(ID_class3)))

write.table(GSEA, file = 'data.txt', append = FALSE, 
            quote = FALSE, sep = '\t',
            row.names = FALSE,
            col.names = TRUE 
)

write(cls, file = 'cls.cls', append = FALSE, 
      #quote = FALSE, 
      sep = '\t',
      #row.names = FALSE,
      #col.names = TRUE, 
      ncolumns = length(cls)
)
