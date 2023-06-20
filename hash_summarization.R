# Written by CS 2023
#Code to look at and interpret hashing results for experiments that use nuclear oligo hashing

# Adapted from analyze_hash_table.R
suppressPackageStartupMessages({
  library(MASS)
  library(reshape2)
  library(VGAM)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(Matrix)
  library(IRdisplay)
  library(tibble)
  library(ggridges)
  library(monocle3)
})

# wrangle paths
setwd("B:/GitHub/behavioral_genomics/hash_ladder") # where the code is located
data_dir <- "G:\\My Drive\\#Projects\\Project_Meningioma\\Data_directory\\Data_Sequencing-results\\20230522_pilot-run-1\\demultiplexed_output" # where the data is located
hashtable.path = file.path(data_dir,"/hash-final/hashTable.out")

hashTable_id <- read.table(hashtable.path, header=F) %>%
  dplyr::select(V1, V2, V3, V4, V5) %>%
  filter(V5 > 30)

colnames(hashTable_id) = c("Sample", "Cell", "Hash", "Axis", "Count")

UMIs_cell <- read.table(file.path(data_dir,"/hash-final/hashUMIs.per.cell"), header=F)
UMIs_cell <- UMIs_cell[,2:3]
colnames(UMIs_cell) <- c("Cell", "Total_RNA")

hashTable <- hashTable_id %>% left_join(UMIs_cell, by="Cell") %>% na.omit()
head(hashTable)

cellid_order <- as.data.frame(hashTable$Cell)
colnames(cellid_order) <- "Cell"
UMIs_hash <- aggregate(hashTable$Count, 
                       by=list(Hash = hashTable$Cell),
                       FUN=sum)
colnames(UMIs_hash) = c("Cell", "Total_hash")

hashTable <- hashTable %>% left_join(UMIs_hash, by="Cell") %>% na.omit()
head(hashTable)

hashTable_unique = hashTable %>% 
  distinct(Cell, .keep_all = T)
head(hashTable_unique)
dim(hashTable_unique)



