###
#   File name : PCA2.R
#   Author    : Hyunjin Kim
#   Date      : Aug 9, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform PCA on the dataset
#
#   * The "PCA.R" is only for the samples from the first batch
#     This code is for the combined raw counts from the first and the second batches
#
#   Instruction
#               1. Source("PCA2.R")
#               2. Run the function "performPCA2" - specify the input file (raw count) and output directory
#               3. The PCA plots will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PCA2.R/PCA2.R")
#               > performPCA2(rCntPath="./data/raw_counts_combined.rda",
#                             outputDir="./results/PCA/")
###

performPCA2 <- function(rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/data/raw_counts_combined.rda",
                        outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/results2/PCA/") {
  
  ### load library
  if(!require(sva)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("sva")
    library(sva)
  }
  
  
  ### load dataset
  load(rCntPath)
  
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  normalizeRNASEQwithVST <- function(readCount, filtering=TRUE) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filtering == TRUE) {
      ### Remove rubbish rows - this will decrease the number of rows
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    }
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  
  ### A function to perform PCA and save a plot
  pca_plot <- function(normalizedMat, grp, title, filePath) {
    
    ### load library
    if(!require(ggfortify)) {
      install.packages("ggfortify")
      library(ggfortify)
    }
    
    ### PCA
    pca_result <- prcomp(t(normalizedMat))
    pca_group <- data.frame(pca_result$x, group=grp)
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save as png
    png(filename=filePath, width = 1000, height = 800)
    print(ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
            labs(title=title) +
            geom_text(aes(label=colnames(normalizedMat)),hjust=0, vjust=0) +
            scale_color_manual(values = colors) +
            theme_classic(base_size = 16))
    dev.off()
    
  }
  
  
  ### A function to perform 3D PCA and save a plot
  pca_plot_3d <- function(normalizedMat, grp, title, filePath) {
    
    ### load library
    if(!require(scatterplot3d)) {
      install.packages("scatterplot3d")
      library(scatterplot3d)
    }
    
    ### PCA
    pca_result <- prcomp(t(normalizedMat))
    pca_group <- data.frame(pca_result$x, group=grp)
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save as png
    png(filename=filePath, width = 1000, height = 800)
    s3d <- scatterplot3d(pca_group[,1:3], pch = 16,
                         color = colors[as.numeric(pca_group$group)],
                         cex.symbols = 0,
                         main = title)
    text(s3d$xyz.convert(pca_group[, 1:3]), labels = rownames(pca_group),
         cex= 1.5, col = colors[as.numeric(pca_group$group)])
    legend("bottomright", legend = levels(pca_group$group),
           col = colors, pch = 16, cex = 1.5)
    dev.off()
    
  }
  
  
  ### perform PCA
  pca_plot(normalizedMat = normalizeRNASEQwithVST(rawCnt),
           grp = raw_count_sample_info$Phenotype,
           title = "Keratinocyte Differentiation",
           filePath = paste0(outputDir, "pca_keratinocyte_differentiation.png"))
  
  ### perform 3D PCA
  pca_plot_3d(normalizedMat = normalizeRNASEQwithVST(rawCnt),
              grp = raw_count_sample_info$Phenotype,
              title = "Keratinocyte Differentiation",
              filePath = paste0(outputDir, "3d_pca_keratinocyte_differentiation.png"))
  
  
  ### batch effect correction
  normCnt_batch_corrected <- ComBat(dat = as.matrix(normalizeRNASEQwithVST(rawCnt)),
                                    batch = raw_count_sample_info$Batch,
                                    mod = model.matrix(~1, data = raw_count_sample_info))
  
  
  ### perform PCA with BC-counts
  pca_plot(normalizedMat = normCnt_batch_corrected,
           grp = raw_count_sample_info$Phenotype,
           title = "Keratinocyte Differentiation - BC",
           filePath = paste0(outputDir, "pca_keratinocyte_differentiation_bc.png"))
  
  ### perform 3D PCA with BC-counts
  pca_plot_3d(normalizedMat = normCnt_batch_corrected,
              grp = raw_count_sample_info$Phenotype,
              title = "Keratinocyte Differentiation - BC",
              filePath = paste0(outputDir, "3d_pca_keratinocyte_differentiation_bc.png"))
  
}

