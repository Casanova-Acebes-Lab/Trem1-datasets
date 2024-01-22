# This script can be used to combine the AURORA, RAP and GEICAM Matrix
# and metadata

library(DESseq)

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

# Reading expression matrix
matrix.AURORA <- read.table(as.character(paste0(rsdir,
                            '/AURORA/matrix_raw_aurora.tsv')),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

matrix.RAP <- read.table(as.character(paste0(rsdir,
                            '/RAP/matrix_raw_rap.tsv')),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

matrix.GEICAM <- read.table(as.character(paste0(rsdir,
                            '/GEICAM/matrix_UQN_geicam.tsv')),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

# Reading metadata
metadata.AURORA <- read.table(as.character(paste0(rsdir,
                            '/AURORA/metadata_aurora.tsv')),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

metadata.RAP <- read.table(as.character(paste0(rsdir,
                            '/RAP/metadata_rap.tsv')),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

metadata.GEICAM <- read.table(as.character(paste0(rsdir,
                            '/GEICAM/metadata_geicam.tsv')),
sep = "\t",
row.names = 1,
header = T,
check.names = F)


# Merging AURORA & RAP datasets

matrix.AURORA$rownames <- row.names(matrix.AURORA)
matrix.RAP$rownames <- row.names(matrix.RAP)

matrix.AURA <- merge(matrix.AURORA, matrix.RAP, 
                    by='rownames')

metadata.AURA <- rbind(metadata.AURORA, metadata.RAP)


# Multi dimensional scalling 

# Filtering low counts

counts.CPM <- cpm(matrix.AURA)
thresh <- counts.CPM > 1
keep <- rowSums(thresh) >= 20
counts.keep <- matrix.AURA[keep,]

# Design matrix

design  <-  model.matrix(~0+cross2$dataset)

# Objeto de DESseq

dds <- DESeqDataSetFromMatrix(counts.keep,
colData = cross2,
design = design)

dim(dds)



write.table(matrix.AURA,file = "results/Combinated/matrix.AURORA.RAP.tsv",
            sep = "\t",
            col.names = NA)




