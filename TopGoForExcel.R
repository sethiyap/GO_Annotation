######################################################################################################
##### Get TopGO ontologies using "elim" algorithm for differentially expressed genes from DESeq2. ####
######################################################################################################
######################################################################################################
##### 1. Requires GO background Database for your organism in the format #############################
##### (Tab separated, without header), mentioned below                   #############################
#####Gene1	GO:0005739,GO:0008150,GO:0046872                             #############################                
#####Gene2	GO:0004672,GO:0005575,GO:0006468                             #############################
######################################################################################################
##### 2. Directory of files containing DESeq output (with extension: *_DEG.xlsx) #####################
#####Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	Set1_C	Set1_T #####
#####Gene1	66.1	-0.15	0.5	-0.28	0.7	0.8	126.5	100.9  #####	
#####Gene2	253.2	-0.31	0.36	-0.86	0.3	0.6	235.0	217.5  #####	
######################################################################################################
##### 3. DESeq output should contain two sheets                      #################################
##### a. "up" regulated genes, name "up"                             #################################
##### b. "down" regulated genes, name "down"                         #################################
######################################################################################################
##### Run the script as ##############################################################################
##### runTopGoForDESeq("./","Temperature/org_geneid2go.map")##########################################
######################################################################################################


#####Required Packages#######
rm(list =ls())
library(topGO)
library("Rgraphviz")
library(psych)
library("xlsx")
library("xlsxjars")
library("rJava")
library("tools")


#######Run the script as mentioned above #######
runTopGoForDESeq = function(dir,GO_DB){
          
### Read the database file and process it ######          
geneID2GO <- readMappings(file = GO_DB,sep = "\t")
str(head(geneID2GO))
GO2geneID <- inverseList(geneID2GO)
str(head(GO2geneID))
geneNames <- names(geneID2GO)
cat(head(geneNames))


##### List the files in the directory ############

filenames = list.files(path=dir,pattern = "*_DEG.xlsx", full.names=T, recursive=FALSE)
length(filenames)

          for(i in 1:length(filenames)){
          
          print(filenames[i])
          
          ###Read sheet of upregulated genes in the DESeq output and process it ######
          
          excelReadUP = read.xlsx(filenames[i], sheetName = "up")
          
          excelUP=as.matrix(excelReadUP[,1])
         
          ###Match your genes in to the database######
          
          geneList <- factor(as.integer(geneNames %in% excelUP))
          
          names(geneList) <- geneNames
          str(geneList)
          
          ####### GO enrichment analysis using TopGo ####
          GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
          
          allGO = genesInTerm(GOdata)
          
          SAM_ANOTATION = lapply(allGO,function(x) x[x %in% excelUP] )

          resultElim = runTest(GOdata, algorithm = "elim", statistic = "fisher")
          
          
          ###### Gives top 20 enriched GO terms ########## 
          UP_GO = GenTable(GOdata,elim=resultElim, topNodes=20 )
          
          allGenes = SAM_ANOTATION[UP_GO[[1]]]
          
          allGenes[length(allGenes) > 1] = lapply(allGenes[length(allGenes) > 1], paste0, collapse = ",")
          
          allGenes=(unlist(allGenes))  
          
          UP_GO$Genes=allGenes
          head(UP_GO)
          
          ###Read sheet of downregulated genes in the DESeq output and process it ######
          
          excelReadDOWN = read.xlsx(filenames[i], sheetName  ="down")
          
          excelDOWN=as.matrix(excelReadDOWN[,1])
          
          geneList <- factor(as.integer(geneNames %in% excelDOWN))
          
          names(geneList) <- geneNames
          
          str(geneList)
          
          GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
          
          allGO = genesInTerm(GOdata)
          
          SAM_ANOTATION = lapply(allGO,function(x) x[x %in% excelDOWN] )
          
          resultElim = runTest(GOdata, algorithm = "elim", statistic = "fisher")
          
          DOWN_GO = GenTable(GOdata,elim=resultElim, topNodes=20 )
          
          allGenes = SAM_ANOTATION[DOWN_GO[[1]]]
          
          allGenes[length(allGenes) > 1] = lapply(allGenes[length(allGenes) > 1], paste0, collapse = ",")
          
          allGenes=(unlist(allGenes))  
          
          DOWN_GO$Genes=allGenes
          
          ####### Write the output in the excel #########
          
          write.xlsx(UP_GO,sheet="UP_GO",file=paste(basename(file_path_sans_ext(filenames[i])),"GO.xlsx", sep="_"), row.names = F, append=F)
          write.xlsx(DOWN_GO,sheet="DOWN_GO",file=paste(basename(file_path_sans_ext(filenames[i])),"GO.xlsx", sep="_"), row.names = F, append=T)
          }
}
