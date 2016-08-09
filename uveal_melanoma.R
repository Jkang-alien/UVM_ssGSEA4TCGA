#untar('gdac.broadinstitute.org_UVM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_Clinical.Level_1.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM-TP.Aggregate_AnalysisFeatures.Level_4.2015082100.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0.tar.gz', exdir = '.')
library('RTCGAToolbox')
library('Cairo')
library('dplyr')
library(survival)
library(rms)
library(ggplot2)
require(gridExtra)

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
immune_mRNA <- data.frame(ID = colnames(mRNA),
                          TSLP = as.numeric(mRNA[rownames(mRNA) == 'TSLP', ]),
                          TNFSF4 = as.numeric(mRNA[rownames(mRNA) == 'TNFSF4', ]),
                          #TSLPR = as.numeric(mRNA[rownames(mRNA) == 'TSLPR', ]),
                          GATA3 = as.numeric(mRNA[rownames(mRNA) == 'GATA3', ]),
                          TBX21 = as.numeric(mRNA[rownames(mRNA) == 'TBX21', ]),
                          IFNG = as.numeric(mRNA[rownames(mRNA) == 'IFNG', ]),
                          IL1B = as.numeric(mRNA[rownames(mRNA) == 'IL1B', ]),
                          TNF = as.numeric(mRNA[rownames(mRNA) == 'TNF', ]),
                          IL12A = as.numeric(mRNA[rownames(mRNA) == 'IL12A', ]),
                          IL12B = as.numeric(mRNA[rownames(mRNA) == 'IL12B', ]),
                          IL10 = as.numeric(mRNA[rownames(mRNA) == 'IL10', ]),
                          IL6 = as.numeric(mRNA[rownames(mRNA) == 'IL6', ]),
                          IL4 = as.numeric(mRNA[rownames(mRNA) == 'IL4', ]),
                          IL5 = as.numeric(mRNA[rownames(mRNA) == 'IL5', ]),
                          IL13 = as.numeric(mRNA[rownames(mRNA) == 'IL13', ]),
                          CD80 = as.numeric(mRNA[rownames(mRNA) == 'CD80', ]),
                          CCL17 = as.numeric(mRNA[rownames(mRNA) == 'CCL17', ]),
                          CCL22 = as.numeric(mRNA[rownames(mRNA) == 'CCL22', ]))

immune_mRNA$Del3p <- immune_mRNA$ID %in% ID_3pDel

IFNG <- ggplot(immune_mRNA,
       aes(x=Del3p, y=IFNG))+ 
 geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, #dotsize=1.0,
               binwidth = 0.6) +
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Absence", "Presence")) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  scale_fill_grey() +
  theme_classic()+
  labs(title="Interferon gamma",x="Loss of Chromosome3", 
      y = "mRNA expression (RSEM)") +
  geom_path(x=c(1,1,2,2),y=c(25,26,26,25))+
  geom_path(x=c(2,2,3,3),y=c(37,38,38,37))+
  geom_path(x=c(3,3,4,4),y=c(49,50,50,49))+
  annotate("text",x=1.5,y=27,label="p=0.012")+
  annotate("text",x=2.5,y=39,label="p<0.0001")+
  annotate("text",x=3.5,y=51,label="p<0.0001")



IL12A <- ggplot(immune_mRNA, 
              aes(x=Del3p,
                  y=IL12A)) +
       geom_dotplot(binaxis='y', stackdir='center',
                    stackratio=1.5, #dotsize=1.0
                    binwidth = 0.05
                    ) +
         scale_x_discrete(breaks=c(FALSE, TRUE),
                          labels=c("Absence", "Presence")) +
       stat_summary(fun.y=median, geom="point", shape=18,
                    size=2, color="red") +
       scale_fill_grey() +
       theme_classic() +
       labs(title="Interleukin12A",
            x="", 
            y = "")

IL12B <- ggplot(immune_mRNA, 
                 aes(x=Del3p,
                     y=IL12B)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, #dotsize=1.0
               binwidth = 0.05
  ) +
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Absence", "Presence")) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  scale_fill_grey() +
  theme_classic() +
  labs(title="Interleukin12B",
       x="", 
       y = "")

CairoPDF(file = './Figures/Th1.pdf', width = 12, height = 4,
         pointsize=18)
multiplot(IFNG, IL12A, IL12B, cols = 3)

dev.off()


########## Th2 ###########################
TSLP <- ggplot(immune_mRNA,
               aes(x=Del3p, y=TSLP))+ 
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, #dotsize=1.0,
               binwidth = 1.0) +
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Absence", "Presence")) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  scale_fill_grey() +
  theme_classic()+
  labs(title="TSLP",x="Loss of Chromosome3", 
       y = "mRNA expression (RSEM)") 



IL4 <- ggplot(immune_mRNA, 
                aes(x=Del3p,
                    y=IL4)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, #dotsize=1.0
               binwidth = 0.005
  ) +
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Absence", "Presence")) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  scale_fill_grey() +
  theme_classic() +
  labs(title="Interleukin4",
       x="", 
       y = "")

IL13 <- ggplot(immune_mRNA, 
                aes(x=Del3p,
                    y=IL13)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, #dotsize=1.0
               binwidth = 0.01
  ) +
  scale_x_discrete(breaks=c(FALSE, TRUE),
                   labels=c("Absence", "Presence")) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  scale_fill_grey() +
  theme_classic() +
  labs(title="Interleukin13",
       x="", 
       y = "")

CairoPDF(file = './Figures/Th2.pdf', width = 12, height = 4,
         pointsize=18)
multiplot(TSLP, IL4, IL13, cols = 3)

dev.off()

wilcox.test(IFNG ~ Del3p, immune_mRNA)
wilcox.test(IL12A ~ Del3p, immune_mRNA)
wilcox.test(IL12B ~ Del3p, immune_mRNA)
wilcox.test(TSLP ~ Del3p, immune_mRNA)
wilcox.test(IL4 ~ Del3p, immune_mRNA)
wilcox.test(IL13 ~ Del3p, immune_mRNA)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




require(gridExtra)
library(grid)

grid.arrange(IFNG, IL4)

pl <- list(IFNG, IL4)
ml <- marrangeGrob(pl, nrow=1, ncol=2)


multiplot(IFNG, IL4, cols = 2)

CairoPDF(file = './Figures/IL_3pDel_subset.pdf', width = 12, height = 16,
         font = 8)
#par(mfrow = c(4,4))
grid.arrange(plot1, plot2)
dev.off()


CairoPDF(file = './Figures/GATA_TBX21.pdf', width = 6, height = 6,
         font = 10)
par(mfrow = c(1,1))

plot(immune_mRNA$TBX21, immune_mRNA$GATA3, 
     col = immune_mRNA$Del3p + 1,
     xlab = 'TBX21 mRNA expression (RSEM)',
     ylab = 'GATA3 mRNA expression (RSEM)'
     )

#abline (h= 10)
#abline (v = 10)

legend(70, 40, c('3pNor', '3pDel'),
       pch = 1,
       col = c(1,2),
       bty = "0")

legend(70, 20, 'Correlation coefficient: 0.82',
       bty = "n")

#abline(0, 1.2)

dev.off()


CairoPDF(file = './Figures/Th2.pdf', width = 12, height = 6,
         font = 10)
par(mfrow = c(1,2))
plot(immune_mRNA$TSLP, immune_mRNA$GATA3,
     col = immune_mRNA$Del3p + 1,
     xlab = 'TSLP mRNA expression (RSEM)',
     ylab = 'GATA3 mRNA expression (RSEM)')
legend(70, 80, c('3pNor', '3pDel'),
       pch = 1,
       col = c(1,2),
       bty = "0")
legend(70, 60, 'Correlation coefficient: -0.15',
       bty = "n")

plot(immune_mRNA$GATA3, immune_mRNA$TNFSF4,
     col = immune_mRNA$Del3p + 1,
     xlab = 'GATA3 mRNA expression (RSEM)',
     ylab = 'TNFSF4 mRNA expression (RSEM)')
legend(40, 30, c('3pNor', '3pDel'),
       pch = 1,
       col = c(1,2),
       bty = "0")
legend(40, 10, 'Correlation coefficient: 0.84',
       bty = "n")

dev.off()

cor(immune_mRNA$GATA3, immune_mRNA$TSLP)

immune_mRNA$ratio_GATA3_TBX21 <- round(immune_mRNA$GATA3/immune_mRNA$TBX21, digits = 2)
clin <- getData(readData, 'Clinical')
clin$months_surv <- clin$days_to_death
clin$months_surv[is.na(clin$days_to_death)] <- clin$days_to_last_followup[is.na(clin$days_to_death)]
clin$ID <- toupper(rownames(clin))
clin$OS_M <- round(as.numeric(clin$months_surv)/30.4, digits = 1)
#clin$f_ratio_GATA3_TBX21 <- factor(clin$ratio_GATA3_TBX21 > 0.8)

library(survival)
library(rms)


data_clin <- merge(clin, immune_mRNA, by = 'ID', all.y = TRUE)
data_clin$f_ratio_GATA3_TBX21 <- factor(data_clin$ratio_GATA3_TBX21 > 1.2)
data_subset <- data_clin %>%
  filter(Del3p == TRUE ) 

strata <- c('ratio GATA/TBX21 > 1.2',
            'ratio GATA/TBX21 <= 1.2')
library(Cairo)

CairoPDF(file = "./Figures/Survival.pdf",  width = 5, height = 5, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
par(
  #mfrow = c(1,2), 
  mar=c(6,4,2,6), mgp = c(2, 1, 0))


fit = npsurv(Surv(OS_M, vital_status == 1)~data_subset$GATA3 > 40,
               data = data_subset)
fit

diff = survdiff(Surv(OS_M, vital_status == 1)~data_subset$GATA3 > 40,
                data = data_subset)
diff

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:2),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=FALSE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.3, 
         y.n.risk=-0.2, cex.n.risk=0.6, pr=FALSE       
)

legend(10, 0.2, strata, lty = c(1:2), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(10, 0.1, 'P-value = 0.412', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

dev.off()

cox <- coxph(Surv(OS_M, vital_status == 1)~(data_subset$GATA3>10)+(data_subset$TBX21>10),
      data = data_subset)
sink('./Figures/stastics.txt')
coxph(Surv(OS_M, vital_status == 1)~(data_subset$GATA3>10) + (data_subset$TBX21>10),
      data = data_subset)
sink()

coxph(Surv(OS_M, vital_status == 1)~(data_subset$TBX21>10) ,
      data = data_subset)

############# IL expression in 3pDel subset ##################

data_subset_GATA_TBX21 <- data_clin %>%
  #filter(Del3p == TRUE ) %>%
  filter (GATA3 >10) %>%
  filter (TBX21 >10)

library(ggplot2)
require(gridExtra)

CairoPDF(file = './Figures/IL_3pDel_subset.pdf', width = 12, height = 16,
         font = 8)
#par(mfrow = c(4,4))
grid.arrange(plot1, plot2)
for (i in 20:35){
  CairoPDF(file = './Figures/IL_3pDel_subset.pdf', width = 12, height = 16,
           font = 8)
  grid.arrange(plot1, plot2)
print (ggplot(data_subset_GATA_TBX21, 
            aes(x=data_subset_GATA_TBX21$f_ratio_GATA3_TBX21,
                y=data_subset_GATA_TBX21[,i])) + 
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center',
                 stackratio=1.5, dotsize=1.2,
                 binwidth = 3))
  

}
  
dev.off()

CairoPDF(file = './Figures/IL_3pDel_subset.pdf', width = 12, height = 16,
         font = 8)
for (i in 20:35) {
print (ggplot(data_subset_GATA_TBX21, 
          aes(x=f_ratio_GATA3_TBX21,
              y=data_subset_GATA_TBX21[i])) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, dotsize=1.2,
               binwidth = 3) + 
  scale_fill_grey() + 
  theme_classic())

}
dev.off()

summary(lm(IFNG~TBX21 * GATA3 , data_subset))
lm <- lm(IFNG~ TBX21 + (TBX21 %in% GATA3) , data_clin)
summary(lm)
plot(lm)

#mean centering linear regression http://www.statedu.com/term/105360

## GATA3 does not inhibit Th1 activity ################

plot(data_subset$ratio_GATA3_TBX21, data_subset$IFNG)
