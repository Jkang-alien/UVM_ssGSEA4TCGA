#GSE39717
source("https://bioconductor.org/biocLite.R")
biocLite("RTCGAToolbox")
library(Biobase)
library(GEOquery)


data_P1 <- read.csv('Supplementary_Table_S1_Worley_et_al.csv', header = TRUE,
                    na.strings = c('NA','ND'))

## ND = not done, NA = not applicable

summary(data_P1)

sample <- read.csv('sample_GSE39717.csv', header = TRUE)
sample <- subset(sample, Collection == 'Fresh frozen sample from enucleated eye')
sample$ID <- gsub('MM ', '', sample$Title)
sample$ID

clin <- merge(data_P1, sample, by.x = 'MM', by.y = 'ID')
summary(clin)
summary(sample)

ID_Chr3_disomy <- clin$Accession[clin$Chr3_status_aCGH == 'Disomy']
ID_Chr3_monosomy <- clin$Accession[clin$Chr3_status_aCGH == 'Monosomy']
######## https://www.bioconductor.org/packages/3.3/bioc/vignettes/GEOquery/inst/doc/GEOquery.html#series
gse<- getGEO(filename = 'GSE39717_family.soft.gz')


gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform})
head(gsmplatforms)

gsmlist = Filter(function(gsm) {Meta(gsm)$platform=='GPL6098'},GSMList(gse))
length(gsmlist)

Table(gsmlist[[1]])[1:5,]
Columns(gsmlist[[1]])[1:5,]
Table(GPLList(gse)[[1]])[2]

probesets <- Table(GPLList(gse)[[1]])[2]
# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
{tab <- Table(x)
mymatch <- match(probesets,tab$ID_REF)
return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})

require(Biobase)
# go through the necessary steps to make a compliant ExpressionSet
rownames(data.matrix) <- probesets$Search_key
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset2 <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
eset2

#data.matrix <- log2(data.matrix)
data.matrix[1:5,]
dim(data.matrix)
data.matrix <- data.matrix[complete.cases(data.matrix),]
data_GSEA <- c()
data_GSEA <- cbind(data.matrix[, colnames(data.matrix) %in% ID_Chr3_monosomy],
                   data.matrix[, colnames(data.matrix) %in% ID_Chr3_disomy])

sum(colnames(data.matrix) %in% ID_Chr3_monosomy)
sum(colnames(data.matrix) %in% ID_Chr3_disomy)
GSEA <- data.frame(NAME = rownames(data.matrix), 
                   DESCRIPTION = rep('na', dim(data_GSEA)[1]), data_GSEA)

write.table(GSEA, file = 'GSEA_GSE39717.txt', 
            quote = FALSE, row.names = FALSE, sep = '\t')

cls <- c(rep(0, sum(colnames(data.matrix) %in% ID_Chr3_monosomy)),
         rep(1, sum(colnames(data.matrix) %in% ID_Chr3_disomy)))

write(cls, file = 'cls_GSE39717.cls', append = FALSE, 
      #quote = FALSE, 
      sep = '\t',
      #row.names = FALSE,
      #col.names = TRUE, 
      ncolumns = length(cls)
)

