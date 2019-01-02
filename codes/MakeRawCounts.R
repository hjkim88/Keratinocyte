###
#   File name : MakeRawCounts.R
#   Author    : Hyunjin Kim
#   Date      : Aug 7, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make raw count matricies from fastq.gz files
#
#   * This code should be run on Linux OS
#
#   Instruction
#               1. Source("MakeRawCounts.R")
#               2. Run the function "makeRCnt" - specify the input directory (fastq.gz) and output directory
#               3. The raw counts will be generated under the output directory
#
#   Example
#               > source("The_directory_of_MakeRawCounts.R/MakeRawCounts.R")
#               > makeRCnt(fastqgzPath="E:/Joanna/AF1805312_R5/",
#                          referencePath="E:/Reference/hg38.fa",
#                          referenceIdxPath="E:/Reference/hg38.index",
#                          outputDir="./data/")
###

makeRCnt <- function(fastqgzPath="/mnt/e/Joanna/AF1805312_R5/",
                     referencePath="/mnt/e/Reference/hg38.fa",
                     referenceIdxPath="/mnt/e/Reference/hg38.index",
                     outputDir="./data/") {
  
  ### load library
  if(!require(Rsubread)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Rsubread")
    library(Rsubread)
  }
  if(!require(org.Hs.eg.db)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
    library(org.Hs.eg.db)
  }
  
  
  ### build the index of the reference genome
  buildindex(basename=referenceIdxPath, reference=referencePath)
  
  
  ### get a list of fastq files
  fastqFiles <- list.files(fastqgzPath)
  fastqFiles <- fastqFiles[which(endsWith(fastqFiles, ".fastq.gz"))]
  
  
  ### iteratively perform alignment
  for(i in 1:(length(fastqFiles)/2)) {
    align(index=referenceIdxPath,
          readfile1=paste0(fastqgzPath, fastqFiles[2*i-1]),
          readfile2=paste0(fastqgzPath, fastqFiles[2*i]),
          output_file=paste0(fastqgzPath, strsplit(fastqFiles[2*i-1], "_", fixed = TRUE)[[1]][1], ".bam"),
          nthreads = 4)
  }
  
  
  ### get a list of created bam files
  bamFiles <- list.files(fastqgzPath)
  bamFiles <- bamFiles[which(endsWith(bamFiles, ".bam"))]
  
  
  ### iteratively get counts from the bam files
  rawCnt <- NULL
  for(i in 1:length(bamFiles)) {
    counts <- featureCounts(files = paste0(fastqgzPath, bamFiles[i]),
                            annot.inbuilt = "hg38",
                            isPairedEnd = TRUE)
    
    if(is.null(rawCnt)) {
      rawCnt <- counts$counts
    } else {
      rawCnt <- cbind(rawCnt, counts$counts)
    }
  }
  
  
  ### numerize the raw counts (now they are characters)
  rawCnt <- as.data.frame(rawCnt)
  rawCnt[1:ncol(rawCnt)] <- lapply(rawCnt[1:ncol(rawCnt)], function(x) as.numeric(as.character(x)))
  
  
  ### set column names for rawCnt and annotate gene symbols
  colnames(rawCnt) <- substr(bamFiles, 1, nchar(bamFiles)-4)
  map_eg_symbol <- mappedkeys(org.Hs.egSYMBOL)
  list_eg2symbol <- as.list(org.Hs.egSYMBOL[map_eg_symbol])
  rawCnt <- cbind(Entrez_ID=rownames(rawCnt),
                  Gene_Symbol=as.character(list_eg2symbol[rownames(rawCnt)]),
                  rawCnt)
  
  
  ### write out the result
  write.table(rawCnt, file = paste0(outputDir, "raw_counts.txt"),
              sep = "\t", row.names = FALSE)
  save(list = c("rawCnt"), file = paste0(outputDir, "raw_counts.rda"))
  
}

