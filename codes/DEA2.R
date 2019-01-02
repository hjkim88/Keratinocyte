###
#   File name : DEA2.R
#   Author    : Hyunjin Kim
#   Date      : Nov 9, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform DE analysis on the combined dataset
#
#   * The "DEA.R" is only for the samples from the first batch
#     This code is for the combined raw counts from the first and the second batches
#
#   Instruction
#               1. Source("DEA2.R")
#               2. Run the function "performDEA2" - specify the input file (raw counts), comparisons, and output directory
#               3. The DEA result will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DEA2.R/DEA2.R")
#               > performDEA2(rCntPath="./data/raw_counts_combined.rda",
#                             outputDir="./results/DEA/Combined/",
#                             comparisons=list(c("D7", "Control"),
#                                              c("D14", "Control"),
#                                              c("D21", "Control"),
#                                              c("D30", "Control"),
#                                              c("D45", "Control"),
#                                              c("D60", "Control")),
#                             fdrThreshold=1e-10)
###

performDEA2 <- function(rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/data/raw_counts_combined.rda",
                        outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/results2/DEA/",
                        comparisons=list(c("D7", "Control"),
                                         c("D14", "Control"),
                                         c("D21", "Control"),
                                         c("D30", "Control"),
                                         c("D45", "Control"),
                                         c("D60", "Control")),
                        fdrThreshold=1e-10) {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  if(!require(org.Hs.eg.db)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
    library(org.Hs.eg.db)
  }
  if(!require(gplots)) {
    install.packages("gplots")
    library(gplots)
  }
  
  
  ### load dataset
  load(rCntPath)
  
  
  ### A function to perform repetitive DE analysis with DESeq2
  deseqWithComparisons <- function(rCnt, grp, exp_class, ctrl_class, bat_eff=NULL) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
    }
    
    ### make a design matrix for DE analysis
    sampleType <- as.character(grp)
    
    if(is.null(bat_eff)) {
      Coldata <- data.frame(sampleType)
    } else {
      batch_eff <- as.character(bat_eff)
      Coldata <- data.frame(sampleType, batch_eff)
    }
    
    rownames(Coldata) <- colnames(rCnt)
    Coldata$sampleType <- relevel(Coldata$sampleType, ref = ctrl_class)
    
    ### data preparation for DE analysis
    if(is.null(bat_eff)) {
      deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType)
    } else {
      deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType+batch_eff)
    }
    
    deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    
    ### run DE analysis
    dea <- DESeq(deSeqData)
    deresults <- results(dea, contrast = c("sampleType", exp_class, ctrl_class))
    if(length(which(is.na(deresults$pvalue))) > 0) {
      deresults$pvalue[which(is.na(deresults$pvalue))] <- 1
    }
    if(length(which(is.na(deresults$padj))) > 0) {
      deresults$padj[which(is.na(deresults$padj))] <- 1
    }
    deresults <- deresults[order(deresults$padj, na.last = TRUE),]
    
    return(deresults)
  }
  
  
  ### A function to print volcano plot of DE analysis with DESeq2 result
  volPlotWithDeseq <- function(deresult, outputFilePath, pvalue=0.05) {
    
    ### load library
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    deresult$padj[which(is.na(deresult$padj))] <- 1
    volcanoData <- as.data.frame(cbind(deresult$log2FoldChange, -log10(deresult$padj), as.character(as.factor(deresult$padj < pvalue))))
    colnames(volcanoData) <- c("logFC", "logFDR", "Significance")
    volcanoData$logFC <- as.numeric(as.character(volcanoData$logFC))
    volcanoData$logFDR <- as.numeric(as.character(volcanoData$logFDR))
    
    s <- strsplit(outputFilePath, "/")
    f1 <- s[[1]][length(s[[1]])]
    
    ggplot(data=volcanoData, aes(x=logFC, y=logFDR, colour=Significance)) + xlab("log2_Fold_Change") + ylab("-log10(FDR)") + ggtitle(paste("Significant (adj.p < ", pvalue, " ) DE genes -", sum(volcanoData$Significance == "TRUE"))) + theme_classic(base_size = 16) + geom_point(alpha=0.4)
    ggsave(filename = outputFilePath)
  }
  
  
  ### Entrez ID - Gene symbol mapping info
  map_eg_symbol <- mappedkeys(org.Hs.egSYMBOL)
  list_eg2symbol <- as.list(org.Hs.egSYMBOL[map_eg_symbol])
  
  
  ### iteratively perform DE analysis for each comparison
  sigGeneNum <- NULL
  for(i in 1:length(comparisons)) {
    ### DESeq2
    de <- deseqWithComparisons(rCnt = rawCnt,
                               grp = raw_count_sample_info$Phenotype,
                               exp_class = comparisons[[i]][1],
                               ctrl_class = comparisons[[i]][2],
                               bat_eff = NULL)
    
    ### set file name
    fileName <- paste("DESeq2", comparisons[[i]][1], "vs", comparisons[[i]][2], sep = "_")
    
    ### annotate gene symbols for the shortened Entrez ID list
    de <- cbind(Entrez_ID=rownames(de), Gene_Symbol=as.character(list_eg2symbol[rownames(de)]), data.frame(de))
    
    ### numerize the factor columns
    idx <- sapply(de, is.factor)
    de[idx] <- lapply(de[idx], function(x) as.character(x))
    
    ### write out the DE result
    write.xlsx2(de, file = paste0(outputDir, fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
    
    ### Volcano plot
    volPlotWithDeseq(de, paste0(outputDir, fileName, "_volPlot.png"), pvalue = fdrThreshold)
    
    ### save the number of significant genes for the comparison
    sigGeneNum <- c(sigGeneNum, length(which(de$padj < fdrThreshold)))
    names(sigGeneNum)[i] <- paste(comparisons[[i]][1], "vs", comparisons[[i]][2], sep = "_")
  }
  
  
  ### draw a barplot to compare the number of DE genes
  png(filename = paste0(outputDir, "The_number_of_DE_genes_comparison.png"),
      width = 1200, height = 1000, res = 110)
  xx <- barplot(sigGeneNum,
                main = paste0("The number of DE genes with FDR < ", fdrThreshold),
                col = "skyblue", ylim = c(0, 1.2 * max(sigGeneNum)))
  text(x = xx, y = sigGeneNum, label = sigGeneNum, pos = 3, cex = 1.5, col = "black")
  dev.off()
  
}


