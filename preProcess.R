# This script can be used to combine the AURORA, RAP and GEICAM Matrix
# and metadata

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

matrix.AURORA.RAP <- merge(matrix.AURORA, matrix.RAP, by=row.names)



