##########################################
### GO enrichment by topGO
##########################################

#Library
library(topGO)
library(tidyverse)

#Current data
today <- Sys.Date()
print(today, quate=F)

##read output of eggNOG mapper: Go information of all genes
geneID2GO <- readMappings(file = "GeneInfo/orf2go.map", sep = "\t", IDsep = ",")
#check
str(head(geneID2GO))
#gene name
geneNames <- names(geneID2GO)

for (Type in c("RedType", "GreenType")){
  print(Type)
  
  for (i in formatC(1:6, width=2,flag="0")){

    print(i)
    
    #read query genes
    FinalVar <- read.csv(paste("kmeans/ClusterInfo/Clus",i,"_GO_",Type,".csv", sep=""),header=T) %>% 
      dplyr::select(GeneID)
    
    #resetting gene names corresponding to GO table
    QueryGenes <- str_c(as.character(FinalVar[,1]), ".1")
    
    #read as list
    list <- unlist(as.list(QueryGenes))
    #count query genes number in all genes
    geneList <- factor(as.integer(geneNames %in% list)) 
    names(geneList) <- geneNames
    geneList 
    
    #Construct topGOdata
    GOdata <- new("topGOdata", 
              ontology = "BP", #BP or CC or MF
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO, 
              nodeSize = 10) 
    
    ###Descrition of parameter
    ##annFUN: 
    #gene2GO this function is used when the annotations are provided as a gene-to-GOs mapping.
    ##nodeSize: 
    #an integer larger or equal to 1. 
    #This parameter is used to prune the GO hierarchy from the termswhich have less than nodeSize annotated genes 
    #(after the true path rule is applied)
    
    #Enrichment test
    resultFis <- runTest(GOdata, 
                     algorithm = "classic", #classi or elim or weight01(default)
                     statistic = "fisher") 
    
    #GenTable
    num_nodes=length(GOdata@graph@nodes) 
    
    GenTable <- GenTable(GOdata, 
         Pvalue = resultFis,
         #orderBy = "weight", 
         #ranksOf = "classic", 
         topNodes = num_nodes)
    
    #Multiple test
    GenTable <- bind_cols(GenTable, p.adjust(GenTable$Pvalue, method = "BH"))
    colnames(GenTable)[7] <- "Adj.Pvalue"
    
    GenTable.Sig <- GenTable %>% filter(Adj.Pvalue < 0.05)
    
    write.csv(GenTable, paste("kmeans/topGO/Clus",i,"_",Type,"_GenTable_elim.csv", sep=""), row.names = FALSE)
    write.csv(GenTable.Sig, paste("kmeans/topGO/Clus",i,"_",Type,"_GenTable_elim_FDR005.csv", sep=""), row.names = FALSE)
    
  }
  }

