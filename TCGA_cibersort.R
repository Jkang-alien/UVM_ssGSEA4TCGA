library('RTCGAToolbox')


getFirehoseDatasets()
getFirehoseRunningDates(last = NULL)
readData = getFirehoseData (dataset="UVM", runDate="20151101",forceDownload = TRUE,
                            Clinic=TRUE, Mutation=FALSE, Methylation=FALSE, RNAseq2_Gene_Norm=TRUE)

clin = getData(readData, "Clinical")
clin[1:5,]
mRNA <- getData(readData, 'RNASeq2GeneNorm')
mRNA[1:10, 1:10]

sum(duplicated(gsub('-...-...-....-..', '', colnames(mRNA))))
colnames(mRNA) <- gsub('-...-...-....-..', '', colnames(mRNA))

write.table(mRNA, file = 'TCGA_mixture.txt', append = FALSE, 
            quote = FALSE, sep = '\t',
            row.names = TRUE,
            col.names = TRUE 
)

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
ID_3qDel <- data$ID[data$CN_3q_Del <= -0.2]
ID_3pNor <- data$ID[data$CN_3p_Del > -0.2]
ID_3qNor <- data$ID[data$CN_3q_Del > -0.2]

sum(ID_3qNor == ID_3pNor)
sum(ID_3qDel == ID_3pDel)
                    
cibersort_result <- read.csv('CIBERSORT.TCGA_UM.csv')
summary(cibersort_result)

cibersort_result$Input.Sample <- gsub('-','\\.', cibersort_result$Input.Sample)

data_ciber <- merge(data, cibersort_result, by.x = 'ID', by.y = 'Input.Sample')
Chr3 <- factor(data_ciber$CN_3p_Del <= -0.2,
               levels = c(FALSE, TRUE),
               labels = c('Disomy', 'Monosomy'))

colnames(data_ciber)

library(NMF)
library(Cairo)
CairoPDF(file = 'cibersort_TCGA.pdf',
         width =7.5, height = 7.5, pointsize = 16)
#layout(matrix(c(1,1,2,2), ncol = 2, byrow = TRUE),
#       widths = c(1,1),
#       heights = c(276,95)) 
#########################################################################

ann <- data.frame(Chr3 = Chr3)

ah<- aheatmap(data_ciber[,252:273], 
              hclustfun=function(d) hclust(d, method="complete"),
              #Colv = colv,
              annRow = ann#,
              #annCol = ann_col,
              #annColors = ann_colors_e,
              #cex = 2,
              #labRow = rep('',dim(data_hc_e)[1]),
              #labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              #Colv = colv
              #reorderfun = function(d, w) reorder(d, 10)
)
dev.off()

t.test(data_ciber$Macrophages.M2~Chr3)


############## single sample GSEA ###########################

ssGSEA <- read.delim('TCGA_UM.gct', sep = '\t')
cell_type <- as.character(ssGSEA$Description)
ID <- colnames(ssGSEA)[3:82]

ES <- t(ssGSEA[,3:82])
#ES <- scale(ES, center = TRUE, scale = TRUE)
#ES <- sweep(ES, 1, apply(ES, 1, median, na.rm = TRUE))
colnames(ES) <- cell_type
## heatmap(ES)
Chr3 <- rownames(ES) %in% ID_3qNor
Chr3 <- factor(Chr3, 
               levels = c(TRUE, FALSE),
               labels = c('Disomy', 'Monosomy'))

hc = hclust(dist(ES), 'ward.D2')

group_immune <- factor(cutree(hc, k = 5), levels = 1:5,
                       labels = c('C', 'A', 'B', 'D', 'E'))

ann <- data.frame(Chr3 = Chr3, group = group_immune)

ann$group <- as.character(ann$group)

ann$group[(ann$Chr3 == 'Monosomy') & (ann$group == 'C')] <- 'D'
ann$group <- factor(ann$group)
df_ann <- data.frame(ID = rownames(ann), ann)
df_ES <- data.frame(ID = rownames(ES), ES)

data_ES <- merge(df_ann, df_ES, by = 'ID')

data_ssGEEA_cli <- merge(data, data_ES, by = 'ID')

boxplot(Th2.cells ~ group, data_ES)

par (mar=c(5.1, 5, 1, 1.1))
bp <- boxplot(Th2.cells ~ Chr3, data_ES,
              ylab = 'Enrichment Score',
              lab = '',
              #border = "white",
              frame = FALSE,
              names = rep('',2),
              xaxt='n',
              outpch = NA
              #ylim = c(3,8.5)
              
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:2,
     pos = -8000,
     tck = -0.01,
     labels=c('Disomy',
              "Monosomy"))

stripchart(Th2.cells ~ Chr3, data_ES,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 3, 3)
y <- c(7.9, 8.0, 8.0, 7.9)-0.9
lines (x,y)
x <- c(2, 2, 4, 4)
y <- c(7.9, 8.0, 8.0, 7.9)-0.1
lines (x,y)


p_value <- c()
for (i in 4:28){
  p <- summary(aov(data_ES[,i] ~ data_ES$group))[[1]][['Pr(>F)']][1]
  p_value <- append(p_value, p)
}

hist(log(p_value))

colnames(ES)[log(p_value) < -30]

aov <- aov(Th17.cells ~ group, data_ES)
summary(aov)
library(ConsensusClusterPlus)

results = ConsensusClusterPlus(ES,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                               title='consensus',
                               clusterAlg="hc",
                               innerLinkage = "ward.D2",
                               finalLinkage = "ward.D2",
                               distance="euclidean",
                               plot="pdf")

library(NMF)
library(Cairo)
CairoPDF(file = 'ssGSEA_TCGA.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES, 
              hclustfun=function(d) hclust(d, method="ward.D2"),
              #Colv = colv,
              annRow = ann,
              #annCol = ann_col,
              #annColors = ann_colors_e,
              #cex = 2,
              #labRow = rep('',dim(data_hc_e)[1]),
              #labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              Colv = results[[4]]$consensusTree
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()

library(survival)
library(rms)

surv_months <- pmax(data_ssGEEA_cli$CLI_days_to_death, 
     data_ssGEEA_cli$CLI_days_to_last_followup,
     na.rm = TRUE)/30.4

diff = survdiff(Surv(surv_months, CLI_vital_status == 1)~ group, 
                data = data_ssGEEA_cli)
diff

setwd('/home/jun/XIST_UT/Figures/')
svg(file = "Figure3.svg", pointsize = 10,
    width = 7.5 , height = 4,)
layout(matrix(c(1,2), ncol = 2, byrow = TRUE))
par(mar=c(5,3,1,4), mgp = c(2, 1, 0))

fit = npsurv(Surv(surv_months, CLI_vital_status == 1)~ group, 
             data = data_ssGEEA_cli)

strata = levels(data_ssGEEA_cli$group)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:5),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data_ssGEEA_cli$group,
         col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.052', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

fit = npsurv(Surv(surv_months, CLI_vital_status == 1)~ Chr3, 
             data = data_ssGEEA_cli)

strata = levels(data_ssGEEA_cli$Chr3)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:2),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data_ssGEEA_cli$group,
         col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)


boxplot(data_ssGEEA_cli)

colnames(data_ssGEEA_cli)

cox <- coxph(Surv(surv_months, CLI_vital_status == 1)~ Th1.cells, 
             data = subset(data_ssGEEA_cli, Chr3 == 'Monosomy'))

summary(cox)

g_Th2.cells <- 
fit = npsurv(Surv(surv_months, CLI_vital_status == 1)~ Chr3, 
             data = data_ssGEEA_cli)

strata = levels(data_ssGEEA_cli$Th1.cells)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:2),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data_ssGEEA_cli$group,
         col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

