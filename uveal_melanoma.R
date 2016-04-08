#untar('gdac.broadinstitute.org_UVM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_Clinical.Level_1.2016012800.0.0.tar.gz', exdir = '.')

untar('gdac.broadinstitute.org_UVM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz', exdir = '.')

untar('gdac.broadinstitute.org_UVM-TP.Aggregate_AnalysisFeatures.Level_4.2015082100.0.0.tar.gz', exdir = '.')

untar('gdac.broadinstitute.org_UVM-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz', exdir = '.')
data <- read.delim('./gdac.broadinstitute.org_UVM-TP.Aggregate_AnalysisFeatures.Level_4.2015082100.0.0/UVM-TP.transferedsamplefeatures.txt')

data[1:5,1:5]
data_gdac <- t(data[,-1])
colnames(data_gdac) <- data[,1]
rownames(data_gdac) <- colnames(data)[-1]
data <- data.frame(data_gdac[-dim(data_gdac)[1],])

data[, c(1:4, 15:250)] <- sapply(data[, c(1:4, 15:250)], as.character)
data[, c(1:4, 15:250)] <- sapply(data[, c(1:4, 15:250)], as.numeric)

summary(data)
colnames(data)

boxplot(data$CN_3p_Del ~ as.factor(data$CLUS_CN_cNMF))

########### Cluster CN NMF ##########################################

################# RNA for GSEA ######################################
############http://www.hammerlab.org/################################


rna <- read.delim('./gdac.broadinstitute.org_UVM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/UVM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt')
rna <- as.matrix(rna[-1,])
#row_number_XIST <- grep('XIST', rna$Hybridization.REF)
dim(rna)
rna[1,-1]
length(unlist(strsplit(rna[,1], '\\|'))[seq(2,2*(dim(rna)[1]), 2)])
rownames(rna) <- unlist(strsplit(rna[,1], '\\|'))[seq(2,2*(dim(rna)[1]), 2)]
ID <- gsub('.[0-9A-Z]{3}.[0-9A-Z]{3}.[0-9A-Z]{4}.[0-9A-Z]{2}', '', colnames(rna))[-1]
index_dup <- duplicated(ID)
sum(index_dup)
