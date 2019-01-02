###
#   File name : CombineTwoDatasets.R
#   Author    : Hyunjin Kim
#   Date      : Nov 8, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Combine two raw count datasets from Joanna to one
#
#   Instruction
#               1. Source("CombineTwoDatasets.R")
#               2. Run the function "combineDS" - specify the input files (raw counts) and output file directory
#               3. The combined raw counts will be generated in the output file directory
#
#   Example
#               > source("The_directory_of_CombineTwoDatasets.R/CombineTwoDatasets.R")
#               > combineDS(rCnt1Path="./data/raw_counts.rda",
#                           rCnt2Path="./data/raw_counts2.rda",
#                           outputDir="./data/")
###

combineDS <- function(rCnt1Path="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/data/raw_counts.rda",
                      rCnt2Path="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/data/data/raw_counts2.rda",
                      outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/data/") {
  
  ### load datasets
  load(rCnt1Path)
  rawCnt1 <- rawCnt
  load(rCnt2Path)
  rawCnt2 <- rawCnt
  rm(rawCnt)
  
  
  ### get common genes
  common_genes <- intersect(rownames(rawCnt1), rownames(rawCnt2))
  
  
  ### reorganize with the common genes
  rawCnt1 <- rawCnt1[common_genes,]
  rawCnt2 <- rawCnt2[common_genes,]
  
  
  ### combine the two datasets
  rawCnt <- cbind(rawCnt1[,-c(1,2)], rawCnt2[,-c(1,2)])
  
  
  ### make sample info
  raw_count_sample_info <- c(rep(0, 3), rep(7, 3), rep(45, 3), rep(60, 2),
                             rep(14, 3), rep(21, 2), rep(30, 3), 21)
  
  
  ### reorder the raw count matrix
  rawCnt <- rawCnt[,order(raw_count_sample_info)]
  raw_count_sample_info <- raw_count_sample_info[order(raw_count_sample_info)]
  
  
  ### slightly change the sample labels
  raw_count_sample_info <- paste0("D", raw_count_sample_info)
  raw_count_sample_info[1:3] <- rep("Control", 3)
  
  
  ### batch info
  raw_count_sample_info <- data.frame(Phenotype=raw_count_sample_info,
                                      Batch=c(rep("B1", 6), rep("B2", 14)))
  rownames(raw_count_sample_info) <- colnames(rawCnt)
  
  
  ### save in RDA format
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"rawCnt\" has 28395 genes as rows and 20 samples as columns")
    writeLines("The Control and D7 samples are from the first batch and the others are from the second batch")
    writeLines("The \"raw_count_sample_info\" has treatment & batch information of the samples")
    writeLines("They are basically time series samples - e.g., D7 means treatment after 7 days")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save as RDA
  save(list = c("rawCnt", "raw_count_sample_info", "README"), file = paste0(outputDir, "raw_counts_combined.rda"))
  
  
  ### save in text format
  rawCnt <- data.frame(Entrez_ID=rownames(rawCnt), rawCnt, check.names = FALSE)
  rawCnt[,1] <- sapply(rawCnt[,1], as.character)
  write.table(rawCnt, file = paste0(outputDir, "raw_counts_combined.txt"), sep = "\t", row.names = FALSE)
  
}

