#untar('gdac.broadinstitute.org_UVM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_Clinical.Level_1.2016012800.0.0.tar.gz', exdir = '.')

untar('gdac.broadinstitute.org_UVM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz', exdir = '.')

clinical <- read.delim('./clinical_data_portal/nationwidechildrens.org_clinical_patient_uvm.txt')

summary(clinical)

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
