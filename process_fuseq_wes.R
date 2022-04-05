#!/usr/bin/env Rscript


################################
# command information
################################
args = commandArgs(trailingOnly=TRUE)


##### analyze input: Rscript process_fuseq_wes.R in=$inPath sqlite=$gtfSqlite txanno=$txAnnofile params=$paramsFn fusiondb=$fusiondbFn out=$outputDir
inPath=txFastaFile=gtfSqlite=txAnnofile=paramsFn=fusiondbFn=paralogdbFn=outputDir=NA
#get command information
args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args)
for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="in") inPath=res[2]
  if (res[1]=="sqlite") gtfSqlite=res[2]
  if (res[1]=="fusiondb") fusiondbFn=res[2]
  if (res[1]=="paralogdb") paralogdbFn=res[2]
  if (res[1]=="params") paramsFn=res[2]
  if (res[1]=="out") outputDir=res[2]
}


outPrefix=paste0(outputDir,"/FuSeq_WES")


###############################
# Load Libraries
###############################

library(GenomicFeatures)

##############################
# Load sources
###############################

source("/path/to/FuSeq_functions.R")
source("/path/to/detectJunctionBreaks.R")
source("/path/to/doBiologicalFilter.R")
source("/path/to/integrateFusion.R")
#### scripts for WES
source("/path/to/processSplitRead_WES.R") #split read processing designed for WES
source("/path/to/postProcessSplitRead_WES.R") #split read post-processing designed for WES
source("/path/to/processFEQ_WES.R") # designed for WES
source("/path/to/processMappedRead_WES.R") # designed for WES
source("/path/to/getSeqLenInfo.R") # designed for WES

##############################
# Prepare params
##############################


if (is.na(paramsFn)){
  cat("\n-----")
  cat("\nThere is no params file. Default settings will be used.")
  FuSeq.params=list()
  FuSeq.params$readStrands="ALL"
  FuSeq.params$chromRef=c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
  #FuSeq.params$chromRef=c("1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","X","Y")
  FuSeq.params$onlyProteinCodingGenes=TRUE
  FuSeq.params$maxSharedCount=5e-2
  FuSeq.params$minGeneDist=1e5
  FuSeq.params$minJunctionDist=1e5
  FuSeq.params$maxInvertedFusionCount=0
  FuSeq.params$maxMRfusionFc=2
  FuSeq.params$maxMRfusionNum=2
  FuSeq.params$sgtMRcount=10
  FuSeq.params$minMR=2
  FuSeq.params$minNonDupMR=2
  FuSeq.params$minSR=1
  FuSeq.params$minScore=3
  FuSeq.params$keepRData=FALSE
  FuSeq.params$exportFasta=FALSE
}else{
  paramIn=read.table(paramsFn, sep="=", header=FALSE)
  FuSeq.params=list()
  FuSeq.params$readStrands=as.character(paramIn[which(paramIn[,1]=="readStrands"),2])
  FuSeq.params$chromRef=trimws(unlist(strsplit(as.character(paramIn[which(paramIn[,1]=="chromRef"),2]),",")))
  FuSeq.params$onlyProteinCodingGenes=as.logical(as.character(paramIn[which(paramIn[,1]=="onlyProteinCodingGenes"),2]))
  FuSeq.params$maxSharedCount=as.double(as.character(paramIn[which(paramIn[,1]=="maxSharedCount"),2]))
  FuSeq.params$minGeneDist=as.double(as.character(paramIn[which(paramIn[,1]=="minGeneDist"),2]))
  FuSeq.params$minJunctionDist=as.double(as.character(paramIn[which(paramIn[,1]=="minJunctionDist"),2]))
  FuSeq.params$maxInvertedFusionCount=as.double(as.character(paramIn[which(paramIn[,1]=="maxInvertedFusionCount"),2]))
  FuSeq.params$maxMRfusionFc=as.double(as.character(paramIn[which(paramIn[,1]=="maxMRfusionFc"),2]))
  FuSeq.params$maxMRfusionNum=as.double(as.character(paramIn[which(paramIn[,1]=="maxMRfusionNum"),2]))
  FuSeq.params$sgtMRcount=as.double(as.character(paramIn[which(paramIn[,1]=="sgtMRcount"),2]))
  FuSeq.params$minMR=as.double(as.character(paramIn[which(paramIn[,1]=="minMR"),2]))
  FuSeq.params$minNonDupMR=as.double(as.character(paramIn[which(paramIn[,1]=="minNonDupMR"),2]))
  FuSeq.params$minSR=as.double(as.character(paramIn[which(paramIn[,1]=="minSR"),2]))
  FuSeq.params$minScore=as.double(as.character(paramIn[which(paramIn[,1]=="minScore"),2]))
  FuSeq.params$keepRData=as.logical(as.character(paramIn[which(paramIn[,1]=="keepRData"),2]))
  FuSeq.params$exportFasta=as.logical(as.character(paramIn[which(paramIn[,1]=="exportFasta"),2]))
}

###############################
# Annotation files
###############################

anntxdb <- loadDb(gtfSqlite)

if (!is.na(fusiondbFn)){
  load(fusiondbFn)
  cat("\n Loading Mitelman Fusion Database ...")
  Mitelman_fusiondb$refseq12=paste(Mitelman_fusiondb$ucsc5,Mitelman_fusiondb$ucsc3,sep="--")
  Mitelman_fusiondb$refseq21=paste(Mitelman_fusiondb$ucsc3,Mitelman_fusiondb$ucsc5,sep="--")
  fglist=Mitelman_fusiondb$refseq12
}


##########################################
# Processing Fusion equivalence classes
##########################################

cat("\n Processing Mapped Reads ");
FuSeq.MR=processMappedRead(inPath,geneAnno=geneAnno,  anntxdb=anntxdb, geeqMap=NULL,FuSeq.params=FuSeq.params)

cat("\n Processing Split Reads ");
FuSeq.SR=processSplitRead(inPath,geneAnno=geneAnno, anntxdb=anntxdb, FuSeq.params=FuSeq.params, txFastaFile=NULL,brolThres=5)


##########################################
# Filters on MR Fusion Candidates
##########################################

# 1. Remove the fusions with shared count
mr=FuSeq.MR$myFusionFinal
pick=which(mr$supportCount==mr$correctedCount)
mr=mr[pick,]

# 2. Remove the ones that two genes are too close
pick=which(mr$chrom1==mr$chrom2 & mr$geneDist<FuSeq.params$minGeneDist)
mr=mr[-pick,]

write.table(mr, file=paste(outPrefix, "_MR_fge.txt", sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)

if (!is.na(fusiondbFn)){
  # 3. Keep only the ones belong to the fusion database
  pick=mr$fusionName %in% Mitelman_fusiondb$refseq12 | mr$fusionName %in% Mitelman_fusiondb$refseq21
  mr_fdb=mr[pick,]
  write.table(mr_fdb, file=paste(outPrefix, "_MR_fge_fdb.txt", sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)
}

#########################################
# Filters on SR Fusion Candidates
#########################################

### Post processing Split Reads
FuSeq.SR.postPro=postProcessSplitRead(inPath, anntxdb, FuSeq.SR, FuSeq.MR, txFastaFile=NULL, FuSeq.params)

sr=FuSeq.SR$myFusionFinal
sr=FuSeq.SR.postPro$myFusionFinal
sr$supportCount=as.integer(sr$supportCount)

### Filter process
myFusionFinal=sr

# 1. At least x read count
myFusionFinal=myFusionFinal[myFusionFinal$supportCount>=1,]
#dim(myFusionFinal)

###
pick=myFusionFinal$tx12Count == myFusionFinal$adjtx12Count
sum(pick)
myFusionFinal=myFusionFinal[pick,]


#keep single results
fgeSet=unique(myFusionFinal$name12)
keepID=NULL
for (i in 1:length(fgeSet)){
  myID=which(myFusionFinal$name12==fgeSet[i])
  keepID=c(keepID,myID[which.max(myFusionFinal$totalCount[myID])])
}
myFusionFinal=myFusionFinal[keepID,]
#dim(myFusionFinal)

write.table(sr, file = paste(outPrefix, "_SR_fge.txt", sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)

if (!is.na(fusiondbFn)){
  pick=myFusionFinal$fusionName %in% Mitelman_fusiondb$refseq12 | myFusionFinal$fusionName %in% Mitelman_fusiondb$refseq21
  sr_fdb=myFusionFinal[pick,]
  #dim(sr_fdb)
  write.table(sr_fdb, file = paste(outPrefix, "_SR_fge_fdb.txt", sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)
}




fusion.MR <- mr
fusion.SR <- sr

fusion.MR$supportCount <- fusion.MR$supportCount + fusion.MR$name21Count

selected_columns = c("fusionName", "gene1", "gene2", "name21", "chrom1", "strand1", "chrom2", "strand2", "supportCount")

myFusionFinal = fusion.SR[,selected_columns]
if (length(fusion.MR$supportCount)>0){
  newID=fusion.MR$fusionName %in% myFusionFinal$fusionName
  myFusionFinal=rbind(myFusionFinal,fusion.MR[!newID,selected_columns])
}

myFusionFinal$MR=myFusionFinal$SR=0
matchID=match(as.character(myFusionFinal$fusionName),as.character(fusion.SR$fusionName))
pick=which(!is.na(matchID))
myFusionFinal$SR[pick]=fusion.SR$supportCount[matchID[pick]]
matchID=match(as.character(myFusionFinal$fusionName),as.character(fusion.MR$fusionName))
pick=which(!is.na(matchID))
myFusionFinal$MR[pick]=fusion.MR$supportCount[matchID[pick]]

myFusionFinal$supportCount=myFusionFinal$SR + myFusionFinal$MR

if (!is.na(fusiondbFn)){
  pick=myFusionFinal$fusionName %in% Mitelman_fusiondb$refseq12 | myFusionFinal$fusionName %in% Mitelman_fusiondb$refseq21
  myFusionFinal$fusionDB=0
  myFusionFinal$fusionDB[pick]=1
}

######################## paralogs

if (!is.na(paralogdbFn)){
  cat("\n Checking paralogs ...")
  isParalogs <- function(gene1, gene2){
    gene1_paralogs <- paralogs[which(gene1 == paralogs$Gene.name), "Human.paralogue.associated.gene.name"]
    return(any(grepl(gene2, gene1_paralogs)))
  }
  load(paralogdbFn)
  myFusionFinal$isParalog <- mapply(isParalogs, myFusionFinal$gene1, myFusionFinal$gene2)
}

cat("\n Exporting results to file  ...")
write.table(myFusionFinal, file = paste(outPrefix, "_FusionFinal.txt", sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)

cat("\n Done!")
