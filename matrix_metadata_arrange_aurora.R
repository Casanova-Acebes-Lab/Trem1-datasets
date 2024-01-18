# Scrip for pre-processing of matrix and metadata from
# AURORA, RAP and GEICAM datasets.

#loading data
matrix_aurora <- read.table(file="input/AURORA/GSE209998_AUR_129_raw_counts.txt",
                            sep = "\t",
                            header = T,
                            check.names = F)

matrix_rap <- read.table(file="input/RAP/GSE193103_salmon_gene.matrix_RAP101_plus_Normals24.txt",
                         sep = "\t",
                         header = T,
                         check.names = F)

matrix_GEICAM <- read.table(file="input/GEICAM/GSE147322_166_AP206_UQN.final.tsv",
                            sep = "\t",
                            header = T,
                            check.names = F)

#removing transcrip variety info function

matrix_arrange<-function(matrix){
  
  matrix[,1] <- gsub('\\..+$', '', matrix[,1])
  
  #deduplicating genes by CV
  vari <- apply(matrix[2: ncol(matrix)],1, var)
  
  matrix_dedup <- matrix[order(vari), ]
  matrix_dedup <- matrix_dedup[!duplicated(matrix_dedup[,1]), ]
  
  row.names(matrix_dedup)<-matrix_dedup[,1]
  matrix_dedup<-matrix_dedup[,-1]
  
  return(matrix_dedup)
}

matrix_aurora <- matrix_arrange(matrix_aurora)
matrix_rap <- matrix_arrange(matrix_rap)
matrix_GEICAM <- matrix_arrange(matrix_GEICAM)

#Arraging Geicam dataset


matrix_GEICAM$gene <- row.names(matrix_GEICAM)
matrix_GEICAM_splt1 <- matrix_GEICAM[1:36,]
matrix_GEICAM_splt1$gene <- row.names(matrix_GEICAM_splt1)
matrix_GEICAM_splt1$gene <- gsub('[|]', '', matrix_GEICAM_splt1$gene)
row.names(matrix_GEICAM_splt1) <- matrix_GEICAM_splt1$gene
matrix_GEICAM_splt1 <- subset(matrix_GEICAM_splt1, select = -c(gene))

matrix_GEICAM_splt2 <- matrix_GEICAM[37:nrow(matrix_GEICAM),]
matrix_GEICAM_splt2$gene <- sub("\\|.*", "", matrix_GEICAM_splt2$gene)

#deduplicating genes by CV in Geicam
matrix <- matrix_GEICAM_splt2 

vari <- apply(matrix[,1:166],1, var)

matrix_dedup <- matrix[order(vari), ]
matrix_dedup <- matrix_dedup[!duplicated(matrix_dedup$gene), ]
matrix_GEICAM_splt2 <- matrix_dedup 

row.names(matrix_GEICAM_splt2) <- matrix_GEICAM_splt2$gene

matrix_GEICAM_splt2 <- subset(matrix_GEICAM_splt2, select = -c(gene))


matrix_geicam <- matrix_GEICAM_splt2

#Saving results for matrix

write.table(matrix_aurora,file = "results/matrix_raw_aurora.tsv",
            sep = "\t",
            col.names = NA)

write.table(matrix_rap,file = "results/matrix_raw_rap.tsv",
            sep = "\t",
            col.names = NA)

write.table(matrix_geicam,file = "results/matrix_UQN_geicam.tsv",
            sep = "\t",
            col.names = NA)



#Creating metadata for Aurora

metadata_aurora <- read.table(file="input/AURORA/metadata_aurora2.tsv",
                              sep = "\t",
                              header = T,
                              check.names = F)

sample <- colnames(metadata_aurora)
tumor <- as.character(metadata_aurora[1,])
tissue <- as.character(metadata_aurora[2,])
treatment <- as.character(metadata_aurora[3,])
time <- as.character(metadata_aurora[4,])
type <- as.character(metadata_aurora[5,])

metadata_aurora <- data.frame(sample)
metadata_aurora$tumor <- tumor
metadata_aurora$tissue <- tissue
metadata_aurora$treatment <- treatment
metadata_aurora$time <- time
metadata_aurora$type <- type



metadata_aurora[,1] <- gsub('human breast cancer total RNA expression', '', 
                            metadata_aurora[,1])
metadata_aurora[,1] <- gsub('human normal brain total RNA expression', '', 
                            metadata_aurora[,1])	
metadata_aurora[,1] <- gsub('human normal brain 2 total RNA expression', '', 
                            metadata_aurora[,1])	
metadata_aurora[,1] <- gsub('human normal lung total RNA expression', '', 
                            metadata_aurora[,1])
metadata_aurora[,1] <- gsub('human normal liver total RNA expression', '', 
                            metadata_aurora[,1])
metadata_aurora[,1] <- gsub('human normal liver total RNA expression', '', 
                            metadata_aurora[,1])
metadata_aurora[,1] <- gsub('human normal breast total RNA expression', '', 
                            metadata_aurora[,1])
metadata_aurora[,1] <- gsub('human normal breast 2 total RNA expression', '', 
                            metadata_aurora[,1])

metadata_aurora[,1] <- gsub("\\[|\\]", "", metadata_aurora[,1])
metadata_aurora$tissue <- gsub("tissue:", "", metadata_aurora$tissue)
metadata_aurora$treatment <- gsub("treatment:", "", metadata_aurora$treatment)
metadata_aurora$time <- gsub("time:", "", metadata_aurora$time)

metadata_aurora$seq <- 'Illumina_HS_2500' 
metadata_aurora$dataset <- 'Aurora'
metadata_aurora$aligner <- 'STAR'
metadata_aurora$quantification <- 'Salmon'

metadata_aurora <- subset(metadata_aurora, select = -c(treatment, time))

metadata_aurora <- metadata_aurora[,-7:-8]

write.table(metadata_aurora,file = "results/metadata_aurora.tsv",
            sep = "\t",
            col.names = NA)

# Creating Metadata for RAP

metadata_rap1 <- read.table(file="input/RAP/metadata_GLP11154.tsv",
                            sep = "\t",
                            header = T,
                            check.names = F)

metadata_rap2 <- read.table(file="input/RAP/metadata_GPL16791.tsv",
                            sep = "\t",
                            header = T,
                            check.names = F)

metadata_rap_arrange <- function(metadata){
  
  sample <- colnames(metadata)
  tumor <- as.character(metadata[2,])
  tissue <- as.character(metadata[1,])
  type <- as.character(metadata[3,])
  seq <- as.character(metadata[7,])
  aligner <- as.character(metadata[5,])
  quantification <- as.character(metadata[6,])
  
  
  metadata <- data.frame(sample)
  metadata$tumor <- tumor
  metadata$tissue <- tissue
  metadata$type <- type
  metadata$seq <- seq
  metadata$aligner <- aligner
  metadata$quantification <- quantification
  
  metadata[,2] <- gsub('cell type:', '', metadata[,2])
  metadata[,2] <- gsub('Metastatic breast cancer', 'Metastatic tumor', 
                       metadata[,2])
  metadata[,2] <- gsub('Primary breast cancer', 'Primary tumor', metadata[,2])
  metadata[,2] <- gsub('Healty tissue', 'Normal tissue', metadata[,2])
  
  metadata[,3] <- gsub('metastasis', '', metadata[,3])
  metadata[,3] <- gsub('Normal', '', metadata[,3])
  metadata[,3] <- gsub('tissue', '', metadata[,3])
  metadata[,3] <- gsub('Primary breast tumor', 'Breast', metadata[,3])
  
  metadata[,4] <- gsub('tissue type:', '', metadata[,4])
  metadata[,4] <- gsub('Formalin-Fixed Paraffin-Embedded', 'FFPE', metadata[,4])
  
  metadata[,6] <- gsub('Maps the data with', '', metadata[,6])
  metadata[,6] <- gsub('STAR 2.7.6a', 'STAR', metadata[,6])
  
  metadata[,7] <- gsub('Quantifies the expression with', '', metadata[,7])
  metadata[,7] <- gsub('Salmon 1.4.0', 'Salmon', metadata[,7])
  
  return(metadata)
}

metadata_rap1 <- metadata_rap_arrange(metadata_rap1)
metadata_rap2 <- metadata_rap_arrange(metadata_rap2)

metadata_RAP <- rbind(metadata_rap1, metadata_rap2)
metadata_RAP$dataset <- 'RAP'

metadata_RAP <- metadata_RAP[,-6:-7]


write.table(metadata_RAP,file = "results/metadata_rap.tsv",
            sep = "\t",
            col.names = NA)

# Building the GEICAM metadata

metadata_geicam <- read.csv(file="input/GEICAM/SraRunTable_geicam.txt",
                            sep = ",",
                            header = T,
                            check.names = F)

metadata_geicam <- data.frame(metadata_geicam$`Library Name`, 
                              metadata_geicam$histological_type, 
                              metadata_geicam$body_site,
                              metadata_geicam$Instrument)
metadata_geicam$type <- 'Fresh frozen'
metadata_geicam$dataset <- 'Geicam'

metadata_geicam$metadata_geicam.body_site <- gsub('Ovarian', 'Ovary', 
                                                  metadata_geicam$metadata_geicam.body_site)
metadata_geicam$metadata_geicam.body_site <- gsub('Node', 'Lymph node', 
                                                  metadata_geicam$metadata_geicam.body_site)
metadata_geicam$metadata_geicam.body_site <- gsub('pericardial effusion', 'Pericardium',
                                                  metadata_geicam$metadata_geicam.body_site)
metadata_geicam$metadata_geicam.body_site <- gsub('Pleural', 'Pleura', 
                                                  metadata_geicam$metadata_geicam.body_site)

colnames(metadata_geicam) <- c('sample', 
                               'tumor', 
                               'tissue', 
                               'seq', 
                               'type', 
                               'dataset')


write.table(metadata_geicam,file = "results/metadata_geicam.tsv",
            sep = "\t",
            col.names = NA)



