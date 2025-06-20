if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("tidyverse",quietly=TRUE))
  BiocManager::install("tidyverse")
if(!require("fgsea",quietly=TRUE))
  BiocManager::install("fgsea")
if(!require("AnnotationDbi",quietly=TRUE))
  BiocManager::install("AnnotationDbi")
if(!require("org.Hs.eg.db",quietly=TRUE))
  BiocManager::install("org.Hs.eg.db")
if(!require("biomaRt",quietly=TRUE))
  BiocManager::install("biomaRt")
library(tidyverse)
library(fgsea)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)

getImportantPaths <- function(can){
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  genes <- getBM(attributes = "entrezgene_id", mart = mart)
  allGenes <- genes[['entrezgene_id']]
  pathways <- gmtPathways("msigdb.v2024.1.Hs.entrez.gmt")
  pathwayNames <- names(pathways)
  A <- str_detect(names(pathways),"MIR")
  B <- str_detect(names(pathways),"SMIRNOV")
  C <- str_detect(names(pathways),"MIRROR")
  D <- str_detect(names(pathways),"LET_7")
  pathwayNamesNoMiR <- names(pathways)[(!A & !D) | B | C]
  davidPaths <- gmtPathways("daivdPathways.gmt")
  biocartaPaths <- gmtPathways("gene_set_library_crisp.gmt")
  KEGGPaths <- gmtPathways("KEGG_hsa.gmt")
  
  allPathways <- c(pathways,davidPaths,biocartaPaths,KEGGPaths)
  allPathNames <- c(pathwayNames,names(davidPaths),names(biocartaPaths),names(KEGGPaths))
  allPathNames <- allPathNames
  print(can)
  allNetGenes <- c()
  allMiRs <- c()
  geneLists <- list()
  miRLists <- list()
  doc <- read.delim(paste(can,"_Modules_Bipartite.txt",sep=""),header=F)
  for(r in rownames(doc)){
    str <- doc[r,]
    vs <- str_split(str,", ")[[1]]
    genes <- vs[!str_detect(vs,"hsa")]
    miRs <- vs[str_detect(vs,"hsa")]
    allNetGenes <- union(allGenes,genes)
    allMiRs <- union(allMiRs,miRs)
    geneLists[[r]] <- genes
    miRLists[[r]] <- miRs
  }
  
  results <- c()
  resultsNoMiR <- list()
  for(r in rownames(doc)){
    print(r)
    modGenes <- as.vector(geneLists[[r]])
    pVals <- c()
    ORs <- c()
    FEs <- c()
    inDF <- c()
    
    pathwayGeneLists <- c()
    
    for(p in allPathNames){
      if(match(p,allPathNames)%%1000 == 0){print(match(p,allPathNames)/length(allPathNames))}
      pathwayGenes <- intersect(allPathways[[p]],allGenes)
      modAndPathwayGenes <- intersect(pathwayGenes,modGenes)
      if(length(modAndPathwayGenes) > 5){
        pathwayGeneLists <- c(pathwayGeneLists,paste(modAndPathwayGenes, collapse=", "))
        a <- length(intersect(modGenes,pathwayGenes))
        b <- length(setdiff(pathwayGenes,modGenes))
        c <- length(setdiff(modGenes,pathwayGenes))
        d <- length(setdiff(allGenes,union(modGenes,pathwayGenes)))
        
        fisher <- fisher.test(rbind(c(a,b),c(c,d)))
        pVals <- c(pVals,fisher$p.value)
        ORs <- c(ORs,fisher$estimate[["odds ratio"]])
        FEs <- c(FEs,(a*length(allGenes)/length(modGenes)/length(pathwayGenes)))
        inDF <- c(inDF,TRUE)
      }
      else{
        inDF <- c(inDF,FALSE)
      }
    }
    padj <- p.adjust(pVals,method="BH")
    fisherDF <- data.frame(allPathNames[inDF],ORs,FEs,pVals,padj,pathwayGeneLists)
    write.csv(fisherDF,paste(can,"_Module",r,"_pathwayAll_FisherResults.csv",sep=""))
    if(nrow(fisherDF) == 0){
      topPaths <- ""
    }
    else{
      topPaths <- (fisherDF %>% filter(FEs > 1, padj < 0.05))
      topPaths <- topPaths[order(topPaths$FEs,decreasing = T),]
      topPaths <- unlist(lapply(rownames(topPaths),function(x){paste(topPaths[x,]$allPathNames.inDF.,"-FE:",topPaths[x,]$FEs,sep = "")}))
    }
    #results[[paste("Module",r,sep="")]] <- topPaths
    results <- c(results,paste(c(paste("Module",r,sep=""),topPaths),collapse=","))
  }
  #resultTab <- paste(results,collapse='\n')
  writeLines(results,sep='\n',con=paste(can,"ModuleBipartitePathwaysAll.csv",sep=""))
}