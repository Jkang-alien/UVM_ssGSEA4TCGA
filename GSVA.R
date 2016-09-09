library(ssGSEA4TCGA)
library(ConsensusClusterPlus)
library(NMF)
library(Cairo)
library(survival)
library(rms)

rdate <- getFirehoseRunningDates(last = NULL)

dset <- getFirehoseDatasets()
dset <- 'UVM'
gs <- gs_gmt('custom.bindea_correct.gmt')

data <- pancancer_ssGSEA(dset, rdate[1], gs)
data_t <- scale(t(data[[1]]), center = TRUE, scale = TRUE)

ann <- data.frame(class = as.factor(results_row[[6]]$consensusClass),
                  type = data[[2]], 
                  Tumor = TvsN)

results_col = ConsensusClusterPlus(data_t,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_col',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

results_row = ConsensusClusterPlus(data[[1]],maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_by_row',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

ann_col <- data.frame(immune_class = factor(results_col[[3]]$consensusClass))
ann <- data.frame(class = as.factor( results_row[[3]]$consensusClass))
rownames(ann) <- gsub('-', '\\.', gsub('-...-...-....-..','',rownames(ann)))

data <- read.delim('./gdac.broadinstitute.org_UVM-TP.Aggregate_AnalysisFeatures.Level_4.2016012800.0.0/UVM-TP.transferedsamplefeatures.txt')

colnames(data)
data[1:5,1:5]
data_gdac <- t(data[,c(-1,-82, -83)])
dim(data_gdac)
colnames(data_gdac) <- data[,1]
rownames(data_gdac) <- colnames(data)[c(-1,-82, -83)]
data <- data.frame(data_gdac)
colnames(data)
data[, c(1:4, 15:251)] <- sapply(data[, c(1:4, 15:251)], as.character)

data[, c(1:4, 15:251)] <- sapply(data[, c(1:4, 15:251)], as.numeric)
data$ID <- rownames(data)

summary(data)
colnames(data)

ID_3pDel <- data$ID[data$CN_3p_Del <= -0.2]
ID_3qDel <- data$ID[data$CN_3q_Del <= -0.2]
ID_3pNor <- data$ID[data$CN_3p_Del > -0.2]
ID_3qNor <- data$ID[data$CN_3q_Del > -0.2]

sum(ID_3qNor == ID_3pNor)
sum(ID_3qDel == ID_3pDel)

Chr3 <- rownames(ann) %in% ID_3qNor
Chr3 <- factor(Chr3, 
               levels = c(TRUE, FALSE),
               labels = c('Disomy', 'Monosomy'))

ann$Chr3 <- Chr3


data <- data.frame(data, ann[match(rownames(ann),data$ID),])
surv_months <- pmax(data$CLI_days_to_death, 
                    data$CLI_days_to_last_followup,
                    na.rm = TRUE)/30.4
data$surv_months <- surv_months

diff = survdiff(Surv(surv_months, CLI_vital_status == 1)~ class, 
                data = data)
diff

diff = survdiff(Surv(surv_months, CLI_vital_status == 1)~ class, 
                data = subset(data, Chr3 == 'Disomy'))
diff

svg(file = "Figure3.svg", pointsize = 10,
    width = 7.5 , height = 4,)
layout(matrix(c(1,2), ncol = 2, byrow = TRUE))
par(mar=c(5,3,1,4), mgp = c(2, 1, 0))

fit = npsurv(Surv(surv_months, CLI_vital_status == 1)~ class, 
             data = subset(data, Chr3 == 'Monosomy'))

strata = levels(data$class)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:5),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data$class,
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
             data = data)

strata = levels(data$Chr3)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:2),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data$class,
         col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)


boxplot(data)

colnames(data)

cox <- coxph(Surv(surv_months, CLI_vital_status == 1)~ Th1.cells, 
             data = subset(data, Chr3 == 'Monosomy'))

summary(cox)
