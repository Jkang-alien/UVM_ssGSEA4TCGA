#GSE39717
source("https://bioconductor.org/biocLite.R")
biocLite("RTCGAToolbox")
library(Biobase)
library(GEOquery)
library(GSVA)
library(ssGSEA4TCGA)
library(ConsensusClusterPlus)
library(NMF)
library(Cairo)
library(survival)
library(rms)



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

######## https://www.bioconductor.org/packages/3.3/bioc/vignettes/GEOquery/inst/doc/GEOquery.html#series
gse<- getGEO(filename = 'GSE39717_family.soft.gz')


gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform})
head(gsmplatforms)

gsmlist = Filter(function(gsm) {Meta(gsm)$platform=='GPL6098'},GSMList(gse))
length(gsmlist)

tab <- Table(gsmlist[[2]])
Columns(gsmlist[[1]])[1:5,]
Table(GPLList(gse)[[1]])[2][1:10,]

probesets <- Table(GPLList(gse)[[1]])

# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs

data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
{tab <- Table(x)
mymatch <- match(GPLLIST_table$ID,tab$ID_REF)
return(tab$VALUE[mymatch])
}))

#data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})

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

gs <- gs_gmt('custom.bindea_correct.gmt')
data <- gsva(data.matrix, gs, method = 'ssgsea')
data_t <- scale(t(data), center = TRUE, scale = TRUE)

pca <- prcomp(data_t,
              center = FALSE,
              scale = FALSE) 

#Chr3[match(clin$Accession, rownames(data_t))] <- clin$Chr.3.status._aCGH
clin$Chr.3.status._aCGH[is.na(clin$Chr.3.status._aCGH)] <- clin$Chr.3.status.FISH[is.na(clin$Chr.3.status._aCGH)]
ID_Chr3_disomy <- clin$Accession[clin$Chr.3.status._aCGH == 'Disomy']
ID_Chr3_monosomy <- clin$Accession[clin$Chr.3.status._aCGH == 'Monosomy']
Chr3 <- rep(NA, 41)
Chr3[rownames(data_t) %in% ID_Chr3_disomy] <- 'Disomy'
Chr3[rownames(data_t) %in% ID_Chr3_monosomy] <- 'Monosomy'


library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(pca,
              obs.scale = 1, var.scale = 1, 
              groups = Chr3, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

CairoPDF(file = 'GSE39717_PCA.pdf',
         width =7.5, height = 7.5, pointsize = 16)
print(g)
dev.off()

results_col = ConsensusClusterPlus(data_t,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_col_gse9717',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

results_row = ConsensusClusterPlus(data,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_row_gse9717',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

ann_col <- data.frame(immune_class = factor(results_col[[3]]$consensusClass))
ann <- data.frame(class = as.factor( results_row[[3]]$consensusClass),
                  Chr3 = Chr3)

CairoPDF(file = 'GSE39717_clustering.pdf',
         width =7.5, height = 7.5, pointsize = 16)
aheatmap(data_t,
         hclustfun=function(d) hclust(dist(d, method = 'euclidean'), method = "ward.D2"),
         annRow = ann,
         annCol = ann_col,
         Colv = results_col[[3]]$consensusTree,
         Rowv = results_row[[3]]$consensusTree,
         labRow = rep('',dim(data_t)[1]))
dev.off()
              