###
#   File name : PathwayAnalysis2.R
#   Author    : Hyunjin Kim
#   Date      : Nov 9, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform pathway analysis on the combined dataset
#
#   * The "PathwayAnalysis.R" is only for the samples from the first batch
#     This code is for the combined raw counts from the first and the second batches
#
#   Instruction
#               1. Source("PathwayAnalysis2.R")
#               2. Run the function "performPA2" - specify the inputs and output directory
#               3. The pathway result will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PathwayAnalysis2.R/PathwayAnalysis2.R")
#               > performPA2(deDir="./results/DEA/Combined/",
#                            outputDir="./results/Pathway/Combined/",
#                            topGeneNum=500,
#                            regulated=c("all", "up", "down"))
###

performPA2 <- function(deDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/results2/DEA/",
                       outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Joanna/results2/Pathway/",
                       topGeneNum=500,
                       regulated=c("all", "up", "down")) {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  
  ### get all the DE result file names in the deDir path
  f <- list.files(deDir)
  f <- f[which(endsWith(f, ".xlsx"))]
  
  
  # ******************************************************************************************
  # Pathway Analysis with clusterProfiler package
  # Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
  #        org          = organism that will be used in the analysis ["human" or "mouse"]
  #                       should be either "human" or "mouse"
  #        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
  #        title        = title of the pathway figure [character]
  #        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
  #        displayNum   = the number of pathways that will be displayed [numeric]
  #                       (If there are many significant pathways show the few top pathways)
  #        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
  #        dir          = file directory path of the output pathway figure [character]
  #
  # Output: Pathway analysis results in figure - using KEGG and GO pathways
  #         The x-axis represents the number of DE genes in the pathway
  #         The y-axis represents pathway names
  #         The color of a bar indicates adjusted p-value from the pathway analysis
  #         For Pathview Result, all colored genes are found DE genes in the pathway,
  #         and the color indicates log2(fold change) of the DE gene from DE analysis
  # ******************************************************************************************
  pathwayAnalysis_CP <- function(geneList,
                                 org,
                                 database,
                                 title="Pathway_Results",
                                 pv_threshold=0.05,
                                 displayNum=Inf,
                                 imgPrint=TRUE,
                                 dir="./") {
    
    ### load library
    if(!require(clusterProfiler)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("clusterProfiler")
      library(clusterProfiler)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    if(!require(xlsx)) {
      install.packages("xlsx")
      library(xlsx)
    }
    
    
    ### colect gene list (Entrez IDs)
    geneList <- geneList[which(!is.na(geneList))]
    
    if(!is.null(geneList)) {
      ### make an empty list
      p <- list()
      
      if(database == "KEGG") {
        ### KEGG Pathway
        kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
        
        if(is.null(kegg_enrich)) {
          writeLines("KEGG Result does not exist")
        } else if(imgPrint == TRUE){
          if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
            result <- kegg_enrich@result
            description <- kegg_enrich@result$Description
          } else {
            result <- kegg_enrich@result[1:displayNum,]
            description <- kegg_enrich@result$Description[1:displayNum]
          }
          
          if(nrow(kegg_enrich) > 0) {
            p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
              theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
              scale_x_discrete(limits = rev(description)) +
              guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
              ggtitle(paste0("KEGG ", title))
            
            png(paste0(dir, "kegg_", title, "_CB.png"), width = 2000, height = 1000)
            print(p[[1]])
            dev.off()
          } else {
            writeLines("KEGG Result does not exist")
          }
        }
        
        if(!is.null(kegg_enrich)) {
          return(kegg_enrich@result)  
        }
      } else if(database == "GO") {
        ### GO Pathway
        if(org == "human") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else if(org == "mouse") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else {
          go_enrich <- NULL
          writeLines(paste("Unknown org variable:", org))
        }
        
        if(is.null(go_enrich)) {
          writeLines("GO Result does not exist")
        } else if(imgPrint == TRUE) {
          if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
            result <- go_enrich@result
            description <- go_enrich@result$Description
          } else {
            result <- go_enrich@result[1:displayNum,]
            description <- go_enrich@result$Description[1:displayNum]
          }
          
          if(nrow(go_enrich) > 0) {
            p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
              theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
              scale_x_discrete(limits = rev(description)) +
              guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
              ggtitle(paste0("GO ", title))
            
            png(paste0(dir, "go_", title, "_CB.png"), width = 2000, height = 1000)
            print(p[[2]])
            dev.off()
          } else {
            writeLines("GO Result does not exist")
          }
        }
        
        if(!is.null(go_enrich)) {
          return(go_enrich@result)  
        }
      } else {
        stop("database prameter should be \"GO\" or \"KEGG\"")
      }
    } else {
      writeLines("geneList = NULL")
    }
  }
  
  
  ### iteratively perform pathway analysis for each comparison
  for(i in 1:length(f)) {
    
    ### load DE result
    de <- read.xlsx2(file = paste0(deDir, f[i]), sheetIndex = 1)
    rownames(de) <- de[,1]
    de <- de[,-1]
    
    ### characterize or numerize the factor columns
    de[1] <- lapply(de[1], function(x) as.character(x))
    de[2:ncol(de)] <- lapply(de[2:ncol(de)], function(x) as.numeric(as.character(x)))
    
    ### regulation setting
    if(regulated[1] == "all") {
      fileName <- paste0(substr(f[i], 1, nchar(f[i])-5))
    } else if(regulated[1] == "up") {
      fileName <- paste0(substr(f[i], 1, nchar(f[i])-5), "_Up")
      de <- de[which(de$log2FoldChange > 0),]
    } else if(regulated[1] == "down") {
      fileName <- paste0(substr(f[i], 1, nchar(f[i])-5), "_Down")
      de <- de[which(de$log2FoldChange < 0),]
    } else {
      stop("ERROR: The \"regulated\" parameter should be \"all\" or \"up\" or \"down\".")
    }
    
    ### filter with the cutoff
    de <- de[order(de$padj),]
    de <- de[1:topGeneNum,]
    
    ### pathway analysis
    ### pathway images
    pKegg <- pathwayAnalysis_CP(geneList = rownames(de), org = "human", database = "KEGG",
                                title = fileName,
                                displayNum = 50, imgPrint = TRUE, dir = outputDir)
    pGo <- pathwayAnalysis_CP(geneList = rownames(de), org = "human", database = "GO",
                              title = fileName,
                              displayNum = 50, imgPrint = TRUE, dir = outputDir)
    ### pathway tables
    if(nrow(pKegg) > 0) {
      write.xlsx2(pKegg, file = paste0(outputDir, fileName, "_KEGG_Table.xlsx"),
                  sheetName = paste0(fileName, "_KEGG"))
    }
    if(nrow(pGo) > 0) {
      write.xlsx2(pGo, file = paste0(outputDir, fileName, "_GO_Table.xlsx"),
                  sheetName = paste0(fileName, "_GO"))
    }
    
  }
  
}
