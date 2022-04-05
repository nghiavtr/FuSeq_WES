############################################################
##### process split reads
############################################################

processSplitRead <-function(inPath,geneAnno, anntxdb, txFastaFile, FuSeq.params, brolThres=10){
  cat("\n ------------------------------------------------------------------")
  cat("\n Processing split reads (SR) from dataset: ",inPath, " read strands:", FuSeq.params$readStrands)
  if (FuSeq.params$keepRData){
    if (is.null(FuSeq.params$outputDir)){ 
      FuSeq.params$outputDir=inPath
      cat("\n No output is set, data will be saved to ", inPath)
    }
  }
  
  ##### get fragment information
  # fragmentInfo=read.csv(paste(inPath,"/fragmentInfo.txt",sep=""), header =TRUE, sep="\t")
  # fragDist = read.table(paste(inPath,"/fragmentDist.txt",sep=""))
  # fragRg=fragDist[fragDist[,2]>0,1]
  # flen.min=min(fragRg)
  # flen.max=max(fragRg)
  ##### find split reads
  cat("\n Get split reads ...")
  fusionGene=fsizes=NULL
  frfiles=list.files(inPath,paste("splitReadInfo_*",sep=""))
  for (i in 1:length(frfiles)){
    tmpDat=read.csv(paste(inPath,"/",frfiles[i],sep=""), header =FALSE, sep="\t",stringsAsFactors=FALSE)
    fusionGene=rbind(fusionGene,tmpDat)
    fsizes=c(fsizes,nrow(tmpDat))
  }
  fsizeLadder=cumsum(fsizes)
  colnames(fusionGene)=c("header","readType","direction","front_tx","front_gene","front_hitpos","front_querypos","front_len","back_tx","back_gene","back_hitpos","back_querypos","back_len","matchedGene","matchedDirect","matchedPos")
  # AS1 Negative strand annotations split
  
  comma_splitstr <- function(x) {
    return(unlist(strsplit(as.character(x), ","))[1])
  }
  fusionGene$front_tx = unlist(lapply(fusionGene$front_tx, comma_splitstr))
  fusionGene$front_gene = unlist(lapply(fusionGene$front_gene, comma_splitstr))
  fusionGene$front_hitpos = unlist(lapply(fusionGene$front_hitpos, comma_splitstr))
  
  fusionGene$back_tx = unlist(lapply(fusionGene$back_tx, comma_splitstr))
  fusionGene$back_gene = unlist(lapply(fusionGene$back_gene, comma_splitstr))
  fusionGene$back_hitpos = unlist(lapply(as.character(fusionGene$back_hitpos), comma_splitstr))
  
  fusionGene$matchedGene = unlist(lapply(fusionGene$matchedGene, comma_splitstr))
  fusionGene$matchedPos = unlist(lapply(fusionGene$matchedPos, comma_splitstr))
  
  fusionGene$tx12=paste(fusionGene$front_tx,fusionGene$back_tx,sep="--")
  fusionGene$tx21=paste(fusionGene$back_tx,fusionGene$front_tx,sep="--")
  fusionGene$name12=paste(fusionGene$front_gene,fusionGene$back_gene,sep="--")
  fusionGene$name21=paste(fusionGene$back_gene,fusionGene$front_gene,sep="--")
  
  
  
  ############# starting process
  cat("\n Extract other biological information...")
  allgenes=unique(c(unique(fusionGene$back_gene),unique(fusionGene$front_gene)))
  allres=select(anntxdb, keys=allgenes, columns=c("GENEID","TXCHROM","TXSTRAND"), keytype = "GENEID")
  colnames(allres)=c("GENEID","chrom","strand")

  #add few biological information
  res=allres[match(fusionGene$front_gene,allres$GENEID),]
  colnames(res)=c("GENEID","chrom1","strand1")
  fusionGene=cbind(fusionGene,res[,-1])
  res=allres[match(fusionGene$back_gene,allres$GENEID),]
  colnames(res)=c("GENEID","chrom2","strand2")
  fusionGene=cbind(fusionGene,res[,-1])
  #filter by chromosomes
  fusionGene=chromFilter(fusionGene,chromRef=FuSeq.params$chromRef)
  fusionGene$gene1=fusionGene$front_gene
  fusionGene$gene2=fusionGene$back_gene


  # matchID=match(fusionGene$gene1,geneAnno[,6])
  # res=geneAnno[matchID,]
  # colnames(res)=paste(colnames(res),"1",sep="")
  # fusionGene=cbind(fusionGene,res[,c(2,4)])
  # matchID=match(fusionGene$gene2,geneAnno[,6])
  # res=geneAnno[matchID,]
  # colnames(res)=paste(colnames(res),"2",sep="")
  # fusionGene=cbind(fusionGene,res[,c(2,4)])
  # #filter by protein coding
  # keepID=which(fusionGene$geneType1=="protein_coding" & fusionGene$geneType2=="protein_coding")
  # fusionGene=fusionGene[keepID,]

  
  ##### get supporting count
  res=table(fusionGene$name12)
  fusionGene$supportCount=res[match(fusionGene$name12,names(res))]

  ##### count of name21
  matchID=match(fusionGene$tx12, fusionGene$tx21)
  tx21Count=fusionGene$supportCount[matchID]
  tx21Count[is.na(tx21Count)]=0
  fusionGene$tx21Count=tx21Count
  
  ##### shared count
  res=table(as.character(fusionGene$header))
  length(res)
  fusionGene$readProp=1/res[match(as.character(fusionGene$header),names(res))]
  
  res=table(fusionGene$tx12)
  fusionGene$tx12Count=res[match(fusionGene$tx12,names(res))]
  adjtx12Count=tapply(as.double(fusionGene$readProp),as.character(fusionGene$tx12),sum)
  fusionGene$adjtx12Count=adjtx12Count[match(as.character(fusionGene$tx12),names(adjtx12Count))]
  
  #gene distance
  geneDist=computeGeneDistance(fusionGene,anntxdb)
  fusionGene$geneDist=geneDist
  
  if (FuSeq.params$keepRData){
    cat("\n Saving the full set of SR fusion candidates ...")
    save(fusionGene,file=paste(FuSeq.params$outputDir,"/FuSeq_SR_fusionGene_full.RData",sep=""))
  }
  
  
  myFusion=fusionGene
  
  #filter by gene distance
  rmID=which(as.character(myFusion$chrom1)==as.character(myFusion$chrom2) &  myFusion$geneDist <= FuSeq.params$minGeneDist)
  if (length(rmID)>0) myFusion=myFusion[-rmID,]
  
  #filter by overlapping sequences
  # brolThres = 10
  # myFusion$brol=myFusion$front_querypos+myFusion$front_len-myFusion$back_querypos
  # myFusionFP=myFusion[myFusion$brol>brolThres,]
  # myFusion=myFusion[myFusion$brol<=brolThres,]
  # 
  # FPlist=unique(as.character(myFusionFP$name12))
  # matchID=match(as.character(myFusion$name12),FPlist)
  # myFusion=myFusion[which(is.na(matchID)),]

  res=table(myFusion$name12)
  myFusion$supportCount2=res[match(myFusion$name12,names(res))]
  #transcript level
  res=table(myFusion$tx12)
  myFusion$tx12Count2=res[match(myFusion$tx12,names(res))]
  #remove all ftx not satisfying the sequence overlapping
  keepID=which(myFusion$tx12Count2==myFusion$tx12Count)
  myFusion=myFusion[keepID,]
  
  
  # myFusion$flen=rep(0,nrow(myFusion))
  # fwID=which(myFusion$direction==0 | myFusion$direction==3)
  # rcID=which(myFusion$direction==1 | myFusion$direction==4)
  # #forward
  # myFusion$flen[fwID]=myFusion$back_querypos[fwID] +  as.numeric(myFusion$matchedPos[fwID])- as.numeric(myFusion$back_hitpos[fwID]) + fragmentInfo$readlen
  # #rc
  # myFusion$flen[rcID]=(as.numeric(myFusion$front_hitpos[rcID])-myFusion$front_querypos[rcID]- as.numeric(myFusion$matchedPos[rcID]) ) + fragmentInfo$readlen
  # myFusion=myFusion[myFusion$flen>=flen.min,]
  # myFusion=myFusion[myFusion$flen<=flen.max,]
  # #mySR$front_hitpos+mySR$front_len-mySR$front_querypos
  # mytest=pnorm(myFusion$flen, mean=fragmentInfo$fragLengthMean, sd=fragmentInfo$fragLengthSd)
  # myFusion$flenTest=mytest
  # myFusion=myFusion[myFusion$flenTest>=0.001,]
  # myFusion=myFusion[myFusion$flenTest<=0.999,]

  ### find duplicate tx
  myDup=duplicated(myFusion$tx12)
  myFusionNoDup=myFusion[!myDup,]
  
  duptx1=table(myFusionNoDup$front_tx)
  duptx2=table(myFusionNoDup$back_tx)
  myFusion=cbind(myFusion,as.integer(duptx1[match(as.character(myFusion$front_tx),names(duptx1))]))
  myFusion=cbind(myFusion,as.integer(duptx2[match(as.character(myFusion$back_tx),names(duptx2))]))
  colnames(myFusion)[c(ncol(myFusion)-1,ncol(myFusion))]=c("duptx1","duptx2")
  
 
  ####### check canonical splicing sites
  cat("\n Remaining fusion reads: ",nrow(myFusion))
  cat("\n Extract extra information from splicing sites ...")
  #cat("\n Check canonical splicing sites: not applied for WES")
  #
  #####preparing annotation
#  library(Biostrings)
#  fasta = readDNAStringSet(txFastaFile)
#  fasta_txnames=sapply(names(fasta), function(x) unlist(strsplit(x," "))[1])
  
  shrinkLen=5
  myFusion$front_brpos= as.numeric(myFusion$front_hitpos)+myFusion$front_len-1
  myFusion$brchposEx5=-1
  myFusion$back_brpos=as.numeric(myFusion$back_hitpos)-(myFusion$back_querypos-myFusion$front_querypos-myFusion$front_len)-1
  myFusion$brchposEx3=-1
  
  

  
  
  if (FuSeq.params$keepRData){
    cat("\n Saving SR fusion candidates with extra information...")
    save(myFusion,file=paste(FuSeq.params$outputDir,"/FuSeq_SR_myFusion.RData",sep=""))
  }
  
#  #NOTE: do not apply this filter for WES
#  cat("\n\n Checking possible paralogs...")
#  gp=paste(geneParalog[,1],"-",geneParalog[,2],sep="")  
#  rmID=myFusion$name12 %in% gp | myFusion$name21 %in% gp
#  rmID=!rmID
#  myFusion$paralog=rep(0,nrow(myFusion))
#  myFusion$paralog[rmID]=1
#  myFusion=myFusion[myFusion$paralog==0,]
  

  # res=list(myFusionFinal=myFusion,fusionGene=fusionGene,fragmentInfo=fragmentInfo,fragDist=fragDist, multiTxts=multi_txts)
  res=list(myFusionFinal=myFusion,fusionGene=fusionGene)
	return(res)
}





