#untar('gdac.broadinstitute.org_UVM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_Clinical.Level_1.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM-TP.Aggregate_AnalysisFeatures.Level_4.2015082100.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz', exdir = '.')
library('RTCGAToolbox')
library('Cairo')


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

ID_3pDel <- data$ID[data$CN_3p_Del <= -0.2]
ID_3pNor <- data$ID[data$CN_3p_Del > -0.2]



########### Cluster CN NMF ##########################################

################# RNA for GSEA ######################################
############http://www.hammerlab.org/################################

getFirehoseDatasets()
getFirehoseRunningDates(last = NULL)
readData = getFirehoseData (dataset="UVM", runDate="20151101",forceDownload = TRUE,
                            Clinic=TRUE, Mutation=FALSE, Methylation=FALSE, RNAseq2_Gene_Norm=TRUE)

#clin = getData(readData, "Clinical")
mRNA <- getData(readData, 'RNASeq2GeneNorm')

colnames(mRNA) <- gsub('-...-...-....-..', '', colnames(mRNA))
colnames(mRNA) <- gsub('-', '\\.', colnames(mRNA))
colnames(mRNA)

data_GSEA <- c()
data_GSEA <- cbind(mRNA[, colnames(mRNA) %in% ID_3pDel],
                   mRNA[, colnames(mRNA) %in% ID_3pNor])

GSEA <- data.frame(NAME = rownames(data_GSEA), 
                   DESCRIPTION = rep('na', dim(data_GSEA)[1]), data_GSEA)
cls <- c(rep(0, length(ID_3pDel)),
         rep(1, length(ID_3pNor)))

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

############# mRNA expression comparison #################
log_TSLP <- log(mRNA[rownames(mRNA) == 'TSLP', ])

immune_mRNA <- data.frame(ID = colnames(mRNA),
                          TSLP = as.numeric(mRNA[rownames(mRNA) == 'TSLP', ]),
                          #TSLPR = as.numeric(mRNA[rownames(mRNA) == 'TSLPR', ]),
                          GATA3 = as.numeric(mRNA[rownames(mRNA) == 'GATA3', ]),
                          TBX21 = as.numeric(mRNA[rownames(mRNA) == 'TBX21', ]),
                          IFNG = as.numeric(mRNA[rownames(mRNA) == 'IFNG', ]),
                          IL1B = as.numeric(mRNA[rownames(mRNA) == 'IL1B', ]),
                          TNF = as.numeric(mRNA[rownames(mRNA) == 'TNF', ]),
                          #IL12 = as.numeric(mRNA[rownames(mRNA) == 'IL12', ]),
                          IL10 = as.numeric(mRNA[rownames(mRNA) == 'IL10', ]),
                          IL6 = as.numeric(mRNA[rownames(mRNA) == 'IL6', ]),
                          IL13 = as.numeric(mRNA[rownames(mRNA) == 'IL13', ]),
                          CD80 = as.numeric(mRNA[rownames(mRNA) == 'CD80', ]),
                          CCL17 = as.numeric(mRNA[rownames(mRNA) == 'CCL17', ]),
                          CCL22 = as.numeric(mRNA[rownames(mRNA) == 'CCL22', ]))


CairoPDF(file = './Figures/boxplots', width = 12, height = 9,
         font = 10)
par(mfrow = c(4,3))

for (i in 2:13){
  boxplot(immune_mRNA[,i] ~ immune_mRNA$Del3p,
  main = colnames(immune_mRNA)[i])
} 

dev.off()

boxplot(log(immune_mRNA$IFNG) ~ immune_mRNA$Del3p)
immune_mRNA$Del3p <- immune_mRNA$ID %in% ID_3pDel



