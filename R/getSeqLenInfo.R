############################################################
##### This function is to detect length of sequneces mapped to transcripts, codes are inherited from function detectJunctionBreaks()
getSeqLenInfo <-function(fgeList,inPath,feqInfo,anntxdb, readStrands="ALL", mapLenProp=0.8){
  myFusionFinal=fgeList

  ##### input fusion reads
  #load fusion reads
  cat("\n Get mapped fusion reads")
  frfiles=list.files(inPath,paste(readStrands,"_fusionMappedReadsChunk_*",sep=""))
  fusionRead=NULL;
  fsizes=NULL;
  for (i in 1:length(frfiles)){
    tmpDat=read.csv(paste(inPath,"/",frfiles[i],sep=""), header =FALSE, sep="\t")
    fusionRead=rbind(fusionRead,tmpDat)
    fsizes=c(fsizes,nrow(tmpDat))
  }
  fsizeLadder=cumsum(fsizes)
  colnames(fusionRead)=c("read1","read2","read1Pos","read2Pos","seq1Pos","seq2Pos","seq1Len","seq2Len")
  #dim(fusionRead)
  cat("\n The number of mapped fusion reads: ",nrow(fusionRead))
  #mapping between readID (index in the file) and transcript set (fusionRead$read1)
  fre1=as.character(fusionRead$read1)
  fre1=trimws(fre1)
  fre2=as.character(fusionRead$read2)
  fre2=trimws(fre2)
  fre1=tapply(c(1:length(fre1)),fre1,c)
  fre2=tapply(c(1:length(fre2)),fre2,c)
  

  #load fragment info - this is not neccessary thi moment
  fragmentInfo=read.csv(paste(inPath,"/fragmentInfo.txt",sep=""), header =TRUE, sep="\t")
  readLen=fragmentInfo[1,1]
  #load the feq file
  feqRaw=feqInfo$feqRaw  
  feqRead1=feqRaw[feqRaw$Read==1,]
  feqRead2=feqRaw[feqRaw$Read==2,]
  #get feq-fge map
  feqFgeMap=feqInfo$feqFgeMap
  feq=feqInfo$feq 
  
  cat("\n Preparing other information ...")

  #mapping between feq and transcript set of a single read, so n feq can map to a single tx
  feqFtxMap1=tapply(as.character(feqRead1$Transcript),feqRead1$Feq,c)
  feqRead1Name=sapply(feqFtxMap1,function(x) paste(x,collapse =" "))
  feqFtxMap2=tapply(as.character(feqRead2$Transcript),feqRead2$Feq,c)
  feqRead2Name=sapply(feqFtxMap2,function(x) paste(x,collapse =" "))
  #feqRead1Name and feqRead2Name: names are feq, values are txset: one feq might have more than 1 txset
  #fre1 and fre2: names are txset, values are readID: 1-n


  seq1Len=as.character(fusionRead$seq1Len)
  seq2Len=as.character(fusionRead$seq2Len)
  seq1Len=lapply(seq1Len,function(x) as.integer(unlist(strsplit(x," "))))
  seq2Len=lapply(seq2Len,function(x) as.integer(unlist(strsplit(x," "))))
  #get unique length
  seq1Len=sapply(seq1Len,unique)
  seq2Len=sapply(seq2Len,unique)

  #get some annotations
  exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionFinal$gene1),as.character(myFusionFinal$gene2))), columns=c("TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
  txToGene=select(anntxdb, keys=unique(as.character(feqRaw[,1])), columns=c("GENEID","TXCHROM"), keytype = "TXNAME")

  cat("\n Detect the information of mapped sequences")
  seqLenInfo=list();  
#  for (i in 1:nrow(myFusionFinal)){
  seqLenInfo <- foreach(i=1:nrow(myFusionFinal)) %dopar%{
    #if (i %% 1000 == 0) cat("Rows processed ",i)
    fgeName=as.character(myFusionFinal$name12)[i]
    gene1=as.character(myFusionFinal$gene1[i])
    gene2=as.character(myFusionFinal$gene2[i])
    myfeqID=feqFgeMap[[fgeName]]
    readL.seqLen=readR.seqLen=readID=list();

    #get mreadID
    mreadID=NULL
    for (j in 1:length(myfeqID)){
      lx=feqRead1Name[myfeqID[j]]
      rx=feqRead2Name[myfeqID[j]]
      intersect(fre1[[lx]],fre2[[rx]])
      mreadID=c(mreadID,intersect(fre1[[lx]],fre2[[rx]]))
    }

    #get mapped length
    readL.seqLen=seq1Len[mreadID]
    readR.seqLen=seq2Len[mreadID]
# #old codes
#    for (j in 1:length(myfeqID)){
#      feqName=names(feq[myfeqID[j]])
#      keepID=intersect(fre1[[feqRead1Name[myfeqID[j]]]],fre2[[feqRead2Name[myfeqID[j]]]])
#      #readID[[feqName]]=keepID
#      #get the list of transcripts in the reads
#      xl=unlist(feqFtxMap1[myfeqID[j]])
#      xr=unlist(feqFtxMap2[myfeqID[j]])
#      #find the positions of all transcripts belonging to the genes of the fusion?
#      IDl=which(txToGene$GENEID[match(xl,txToGene$TXNAME)]==gene1)
#      IDr=which(txToGene$GENEID[match(xr,txToGene$TXNAME)]==gene2)
#      #select only one because the others will locate to the same posistion
#      posl=IDl[1]
#      posr=IDr[1]
#      txSeqlenl=sapply(seq1Len[keepID],function(x) x[posl])
#      txSeqlenr=sapply(seq2Len[keepID],function(x) x[posr])
#      readL.seqLen[[feqName]]=txSeqlenl
#      readR.seqLen[[feqName]]=txSeqlenr
#    } #end of j

    res=list(readL.seqLen=readL.seqLen, readR.seqLen=readR.seqLen,mreadID=mreadID)
    #seqLenInfo[[i]]=res;
    res
  }#end of loop

mapLen_thres=readLen*mapLenProp
mapLenCount=sapply(seqLenInfo, function(x){
  x1=unlist(x$readL.seqLen)
  x2=unlist(x$readR.seqLen)
  res=sum(x1>mapLen_thres & x2>mapLen_thres)
  return(res)
})
myFusionFinal$mapLenCount=mapLenCount
myFusionFinal=myFusionFinal[myFusionFinal$mapLenCount>0,]

return(list(myFusionFinal=myFusionFinal,seqLenInfo=seqLenInfo,fsizeLadder=fsizeLadder))
}

