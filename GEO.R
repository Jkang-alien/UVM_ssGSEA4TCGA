########## GEO #####################
source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("GEOquery")
library(Biobase)
library(GEOquery)

sample <- read.csv('sample_subset_GEO.csv')
data_P1 <- read.delim('table_p1.txt', header = FALSE, sep = '')
colnames(data_P1) <- c('Sample_ID', 'Gender', 'Age', 'Thickness',
                       'Size', 'Location', 'Cell type', 'Sdera',
                       'Metastasis', 'DFS_m', 'Chr3', 'Chr8',
                       'PreviousTreat', 'mdda_9/syntenin')

data_P1$ID <- gsub('MU_', '', data_P1$Sample_ID)
sample$ID <- gsub('Uveal melanoma ', '', sample$Title)
data_sample_geo <- merge(sample, data_P1, by = 'ID')
data_P1$ID %in% sample$ID


ID_chr3_m <- data_sample_geo$Samples[data_sample_geo$PreviousTreat == 'none' &
                          data_sample_geo$Chr3 == 'm' ]
ID_chr3_d <- data_sample_geo$Samples[data_sample_geo$PreviousTreat == 'none' &
                                data_sample_geo$Chr3 == 'd' ]

gds4281<- getGEO(filename='GDS4281_full.soft.gz')
Meta(gds4281)$description

data = Table(gds4281)
colnames(data)
data_idenfier <- Table(gds4281)[,c(2,3:31)]

data = Table(gds4281)[,c(1,3:31)]

data_GSEA <- c()
data_GSEA <- cbind(data[, colnames(data) %in% ID_chr3_m],
                   data[, colnames(data) %in% ID_chr3_d])

GSEA <- data.frame(NAME = data$ID_REF, 
                   DESCRIPTION = rep('na', dim(data_GSEA)[1]), data_GSEA)

write.table(GSEA, file = 'GSEA_GEO4281.txt', 
            quote = FALSE, row.names = FALSE, sep = '\t')

cls <- c(rep(0, length(ID_chr3_m)),
         rep(1, length(ID_chr3_d)))

write(cls, file = 'cls_geo4281.cls', append = FALSE, 
      #quote = FALSE, 
      sep = '\t',
      #row.names = FALSE,
      #col.names = TRUE, 
      ncolumns = length(cls)
)

