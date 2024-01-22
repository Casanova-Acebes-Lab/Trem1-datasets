# Scrip for pre-processing of matrix and metadata from
# AURORA, RAP and GEICAM datasets.

library(stringr)

# Reading th config file

config <- read.table("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
rsdir <- values["rsdir"]
datadir <- values["datadir"]

outdir <- paste0(rsdir)

#loading data
matrix_aurora <- read.table(as.character(paste0(datadir,
                            '/AURORA/GSE209998_AUR_129_raw_counts.txt')),
                            sep = "\t",
                            header = T,
                            check.names = F)

matrix_rap <- read.table(as.character(paste0(datadir,
                          "/RAP/GSE193103_salmon_gene.matrix_RAP101_plus_Normals24.txt")),
                          sep = "\t",
                          header = T,
                          check.names = F)

matrix_GEICAM <- read.table(as.character(paste0(datadir,
                            "/GEICAM/GSE147322_166_AP206_UQN.final.tsv")),
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

write.table(matrix_aurora, paste0(outdir,
            "/AURORA/matrix_raw_aurora.tsv"),
            sep = "\t",
            col.names = NA)

write.table(matrix_rap,paste0(outdir,
            "/RAP/matrix_raw_rap.tsv"),
            sep = "\t",
            col.names = NA)

write.table(matrix_geicam,paste0(outdir,
            "/GEICAM/matrix_UQN_geicam.tsv"),
            sep = "\t",
            col.names = NA)



#Creating metadata for Aurora

metadata_aurora <- read.table(as.character(paste0(datadir,
                              "/AURORA/metadata_aurora2.tsv")),
                              sep = "\t",
                              header = T,
                              check.names = F)

# Applying metadata adjustments

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

metadata_aurora$seq <- 'Illumina HiSeq 2500' 
metadata_aurora$dataset <- 'Aurora'
metadata_aurora$aligner <- 'STAR'
metadata_aurora$quantification <- 'Salmon'

metadata_aurora <- subset(metadata_aurora, select = -c(treatment, time))

metadata_aurora <- metadata_aurora[,-7:-8]

metadata_aurora$sample <- str_squish(metadata_aurora$sample)
metadata_aurora$tumor <- str_squish(metadata_aurora$tumor)
metadata_aurora$tissue <- str_squish(metadata_aurora$tissue)
metadata_aurora$type <- str_squish(metadata_aurora$type)
metadata_aurora$seq <- str_squish(metadata_aurora$seq)

write.table(metadata_aurora,paste0(outdir,
            "/AURORA/metadata_aurora.tsv"),
            sep = "\t",
            col.names = NA)

# Creating Metadata for RAP

metadata_rap1 <- read.table(as.character(paste0(datadir,
                            "/RAP/metadata_GLP11154.tsv")),
                            sep = "\t",
                            header = T,
                            check.names = F)

metadata_rap2 <- read.table(as.character(paste0(datadir,
                            "/RAP/metadata_GPL16791.tsv")),
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

metadata_RAP$sample <- str_squish(metadata_RAP$sample)
metadata_RAP$tumor <- str_squish(metadata_RAP$tumor)
metadata_RAP$tissue <- str_squish(metadata_RAP$tissue)
metadata_RAP$type <- str_squish(metadata_RAP$type)
metadata_RAP$seq <- str_squish(metadata_RAP$seq)

write.table(metadata_RAP,paste0(outdir,
            "/RAP/metadata_rap.tsv"),
            sep = "\t",
            col.names = NA)

# Building the GEICAM metadata

metadata_geicam <- read.table(as.character(paste0(datadir, 
                            "/GEICAM/SraRunTable_geicam.txt")),
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
metadata_geicam$metadata_geicam.histological_type <- gsub('Breast cancer Metastasis', 'Metastatic tumor', 
                                                  metadata_geicam$metadata_geicam.histological_type)
metadata_geicam$metadata_geicam.histological_type <- gsub('Breast primary tumor', 'Primary tumor', 
                                                  metadata_geicam$metadata_geicam.histological_type)


colnames(metadata_geicam) <- c('sample', 
                               'tumor', 
                               'tissue', 
                               'seq', 
                               'type', 
                               'dataset')


metadata_geicam$sample <- str_squish(metadata_geicam$sample)
metadata_geicam$tumor <- str_squish(metadata_geicam$tumor)
metadata_geicam$tissue <- str_squish(metadata_geicam$tissue)
metadata_geicam$type <- str_squish(metadata_geicam$type)
metadata_geicam$seq <- str_squish(metadata_geicam$seq)

write.table(metadata_geicam,paste0(outdir,
            "/GEICAM/metadata_geicam.tsv"),
            sep = "\t",
            col.names = NA)


# Merging cohorts

metadata_AU_RA_GE <- rbind(metadata_aurora,
                          metadata_RAP,
                          metadata_geicam)

write.table(metadata_AU_RA_GE, paste0(outdir,
            "/metadata_AU_RA_GE.tsv"),
            sep = "\t",
            col.names = NA)

