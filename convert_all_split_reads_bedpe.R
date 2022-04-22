  frfiles="splitReadInfo.txt"
  fusionGene=read.csv(frfiles, header =FALSE, sep="\t",stringsAsFactors=FALSE)
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
  
  res=table(fusionGene$name12)
  fusionGene$supportCount=unclass(res[match(fusionGene$name12,names(res))])
  
  chrom1=sapply(fusionGene$front_tx,function(x){
    return(unlist(strsplit(as.character(x), "__"))[2])
  },USE.NAMES = FALSE)

  contig1_start=sapply(fusionGene$front_tx,function(x){
    return(as.integer(unlist(strsplit(as.character(x), "__"))[3]))
  },USE.NAMES = FALSE)
  start1=contig1_start+as.integer(fusionGene$front_hitpos)
  end1=start1+as.integer(fusionGene$front_len)

  chrom2=sapply(fusionGene$back_tx,function(x){
    return(unlist(strsplit(as.character(x), "__"))[2])
  },USE.NAMES = FALSE)

 contig2_start=sapply(fusionGene$back_tx,function(x){
    return(as.integer(unlist(strsplit(as.character(x), "__"))[3]))
  },USE.NAMES = FALSE)
  start2=contig2_start+as.integer(fusionGene$back_hitpos)
  end2=start2+as.integer(fusionGene$back_len)
  

  bedpe=data.frame(chrom1=chrom1,start1=start1,end1=end1,chrom2=chrom2,start2=start2,end2=end2,name=fusionGene$name12,score=fusionGene$supportCount, strand1=".", strand2=".")
  bedpe_final=bedpe

  res=unclass(table(bedpe$name))
  minStart1=tapply(bedpe$start1,bedpe$name,min)
  maxEnd1=tapply(bedpe$end1,bedpe$name,max)
  minStart2=tapply(bedpe$start2,bedpe$name,min)
  maxEnd2=tapply(bedpe$end2,bedpe$name,max)
  bedpe_final$start1=minStart1[bedpe$name]
  bedpe_final$end1=maxEnd1[bedpe$name]
  bedpe_final$start2=minStart2[bedpe$name]
  bedpe_final$end2=maxEnd2[bedpe$name]

  bedpe_final=bedpe_final[!duplicated(bedpe_final$name),]


  write.table(bedpe_final, file="splitReadInfo.BEDPE",col.names=TRUE, row.names=FALSE, quote=FALSE)