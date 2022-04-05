##############################
##### extract gtf file for wes exon
## Update: 30Otc2019/Nghia:
# - index fully genes from known fusions in the database if its length less than 500K bases
## Update: 22Otc2019/Nghia:
# - Use RefSeq annotation
# - Prioritise indexing full genes from fusion databases
## Update: 10May2019/Nghia:
# refine the gtf format
## Update: 08Apr2019/Nghia:
# for small genes (<20 000bp, cover 75% genes), we index fully their genome sequences
## Update: 05Apr2019/Nghia:
# Run in turn for individual genes in chromosomes
# For big genes, extend each exon a range of maxDist to both left and right sides, then combind overlapping exons into super contigs
# Divide the super contigs or small genes into smaller contigs with a fixed window size; two consecutive contigs are overlapping in a region of r-1 (99bp)
# Export gtf of the contigs

##############################

args = commandArgs(trailingOnly=TRUE)

#pre-setting
bigGeneLen=130000 # so, 90% genes can be fully indexed
#bigGeneLen=0 # so, all genes are fully indexed
maxDist=1000 #extend 3000bp for each exon
contigSize=3000 #size of divided contigs
bufLen=1000 #merge smaller contigs if it is too small (<bufLen) 
boundaryExtension=1000#extend to boundary of a gene, not exon
readLen=100

##### analyze input: Rscript extract_gtf.R genomefasta=genomeFastaFile sqlite=gtfSqlite fusiondb=fusiondbFn out=gtfoutFn
genomeFastaFile=gtfSqlite=fusiondbFn=gtfoutFn=NA
#get command information
args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args)
for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="genomefasta") genomeFastaFile=res[2]  
  if (res[1]=="sqlite") gtfSqlite=res[2]
  if (res[1]=="fusiondb") fusiondbFn=res[2]
  if (res[1]=="out") gtfoutFn=res[2]
  #other optional parameters
  if (res[1]=="bigGeneLen") bigGeneLen=as.integer(res[2])
  if (res[1]=="maxDist") maxDist=as.integer(res[2])
  if (res[1]=="contigSize") contigSize=as.integer(res[2])
  if (res[1]=="bufLen") bufLen=as.integer(res[2])
  if (res[1]=="boundaryExtension") boundaryExtension=as.integer(res[2])
  if (res[1]=="readLen") readLen=as.integer(res[2])
}


library("GenomicFeatures")
library("Biostrings")



#output file
if (is.na(gtfoutFn)) gtfoutFn=paste("reference_wes_contigSize",contigSize,"_bigLen",bigGeneLen,"_r",readLen,".gtf",sep="")


anntxdb <- loadDb(gtfSqlite)
fasta = readDNAStringSet(genomeFastaFile)
fasta_chrnames_fasta=sapply(names(fasta), function(x) unlist(strsplit(x," "))[1])
chrlen=width(fasta)
names(chrlen)=fasta_chrnames_fasta
#fasta_chrnames_fasta=fasta_chrnames_fasta[1:25] #c(as.character(c(1:22)),"X","Y","MT")
fasta_chrnames_fasta=paste0("chr",c(as.character(c(1:22)),"X","Y","M"))
genes.all=genes(anntxdb,single.strand.genes.only=TRUE) # NOTE: some genes are from different strands in UCSC annotation, but they do not exist in the ensembl annotation
fasta_chrnames=sort(unique(as.character(seqnames(genes.all))))
fasta_chrnames = intersect(fasta_chrnames,fasta_chrnames_fasta)
#fasta_chrnames

load(fusiondbFn)
#dim(Mitelman_fusiondb)


############## select the threshold to define big genes
fgenes=unique(c(as.character(Mitelman_fusiondb$ucsc5),as.character(Mitelman_fusiondb$ucsc3))) #for refseq
x=width(genes.all)
###refine the fgenes
allgn=names(genes.all)
pick=which(allgn %in% fgenes)
allgn_share=allgn[pick]
allgn_len=x[pick]
#fully index only genes of fusion with tx less than 500000
bigFusionGeneLen=max(allgn_len)+1 #consider all

pick=which(allgn_len<bigFusionGeneLen)
allgn_share=allgn_share[pick]
#update fgenes
length(fgenes)
fgenes=fgenes[fgenes %in% allgn_share]
length(fgenes)
sum(allgn_len>=bigFusionGeneLen)

##define a threshold of length for big genes
#bigGeneLen=130000 # so, 90% genes can be fully indexed
##bigGeneLen=0 # so, all genes are fully indexed
######
###ex.wes.fa=DNAStringSet("TCGA") #start with a dumpy DNA sequence
##chrID=1
##maxDist=1000 #extend 3000bp for each exon
##contigSize=3000 #size of divided contigs
##bufLen=1000 #merge smaller contigs if it is too small (<bufLen) 
##boundaryExtension=1000#extend to boundary of a gene, not exon
##readLen=100


#source("getRcodes.R")

for (chrID in 1:length(fasta_chrnames)){
  chr=fasta_chrnames[chrID]
  len=chrlen[names(chrlen)==chr]
  cat("Processing chromosome ",chr)
  #myex=exons(anntxdb,columns=c("EXONID"),filter=list(exon_chrom=chr))
  genes.chr=genes(anntxdb, filter=list(tx_chrom = chr))
  cat(" # of genes ",length(genes.chr))
  genes.exon.map=select(anntxdb, keys=names(genes.chr), columns=c("GENEID","EXONID","EXONSTART","EXONEND","EXONSTRAND","EXONCHROM"), keytype = "GENEID")
  #limit to the current chromosome. NOTE: some genes are annotated in different chromosomes in UCSC hg38
  genes.exon.map=genes.exon.map[genes.exon.map$EXONCHROM %in% chr,]
  #remove exons from different chromosomes: some genes have exons from two chromosomes, e.g. AKAP17A: This requires for genes from UCSC, not Ensembl
  #genes.exon.map=genes.exon.map[genes.exon.map$EXONCHROM==chr,]
  cat(" # of exons ",nrow(genes.exon.map))
  for (g in names(genes.chr)){
#    if((g %in% fgenes)){
    #g=names(genes.chr)[1]
    GEmat=genes.exon.map[genes.exon.map$GENEID==g,]
    GEmat=GEmat[order(GEmat$EXONID),] #sort exons by increasing order for forward strand    
    #get extra boundaryExtension from the boundaries of the gene
    startR=min(GEmat$EXONSTART)-boundaryExtension
    endR=max(GEmat$EXONEND)+boundaryExtension
    startR=ifelse(startR<1,1,startR)
    endR=ifelse(endR>len,len,endR)
    glen= endR-startR +1
    #export gene information to gtf
    exGtf=c(chr, "protein_coding","gene",min(startR), max(endR),".",GEmat$EXONSTRAND[1],".", paste("gene_id ",'"',g,'"; gene_name ','"',g,'"; gene_source "MEB"; gene_biotype "wes";',sep=""))
    exGtf=paste(exGtf,collapse = "\t")
    names(exGtf)=NULL
    write(exGtf, file=gtfoutFn, append = TRUE, sep="\t")
    if (glen > bigGeneLen & !(g %in% fgenes)){ #if it is a big gene, then make supercontigs then divive thems to smaller trunks
#      if (glen < 0){ #not do this anymore
      #get extra maxDist=5000bp
      GEmat$EXONSTART2=GEmat$EXONSTART - maxDist
      GEmat$EXONEND2=GEmat$EXONEND + maxDist
      #re-adjust for the first and last exons
      pick=which(GEmat$EXONSTART %in% startR)
      GEmat$EXONSTART2[pick]=GEmat$EXONSTART2[pick] + maxDist - boundaryExtension
      pick=which(GEmat$EXONEND %in% endR)
      GEmat$EXONEND2[pick]=GEmat$EXONEND2[pick] - maxDist + boundaryExtension         
      #select the coding area inside the region.
      GEmat=GEmat[order(GEmat$EXONSTART2),]
      exstart=GEmat$EXONSTART2/1e6
      exend=GEmat$EXONEND2/1e6  
      #get clusters of exons
      exCluster=rep(-1,length(exstart))
      for (j in 1:length(exstart)){
        if (exCluster[j]==-1){
          exCluster[j]=j
          pick=(exstart[j]-exstart)*(exstart[j]-exend)<=0 | (exend[j]-exstart)*(exend[j]-exend)<=0
          #assign to the cluster with min index
          x=exCluster[pick]
          x=x[x>0]
          x=min(x)
          exCluster[pick]=x
        }
      }
      #update exCluster
      exCluster_u=exCluster
      repeat{
        exCluster=exCluster[exCluster]
        if (sum(exCluster_u!=exCluster)==0) break()
        exCluster_u=exCluster
      }
      #create exon regions from the exon clusters
      clusterID=unique(exCluster)
      startR=endR=NULL #new startR and endR
      for (j in 1:length(clusterID)){
        cID=clusterID[j]
        cEx=GEmat[exCluster==cID,]
        startR=c(startR,min(cEx$EXONSTART2))
        endR=c(endR,max(cEx$EXONEND2))
      }
      startR=ifelse(startR<1,1,startR)
      endR=ifelse(endR>len,len,endR)    
    }#end for big gene

    #divide the super contigs into smaller contigs with a fixed window size
    startSC=startR
    endSC=endR
    startR=endR=NULL
    for (j in 1:length(startSC)){
      SCsize=endSC[j]-startSC[j]+1    
      x=0
      y=SCsize-1
      if (SCsize > (contigSize+bufLen)){
        x=seq(0,SCsize,contigSize)
        y=c(x[-1],SCsize)
        x[-1]=x[-1]-readLen+1
        z=y-x
        t=length(z)
        if(z[t] < bufLen){
          x=x[-t]
          y[t-1]=y[t]
          y=y[-t]
        }
      }
      startR=c(startR,startSC[j]+x)
      endR=c(endR,startSC[j]+y)
    }

    #sort again startR
    myorder=order(startR)
    startR=startR[myorder]
    endR=endR[myorder]
    startR=ifelse(startR<1,1,startR)
    endR=ifelse(endR>len,len,endR)
    startR=as.character(as.integer(startR))
    endR=as.character(as.integer(endR))
    ex.names=paste(g,"__",chr,"__",startR,"__",endR,sep="") #name in gtf
    #export transcript
    exGtf=sapply(c(1:length(startR)),function(t){
      x=c(chr, "processed_transcript","transcript",startR[t], endR[t],".",GEmat$EXONSTRAND[1],".", paste("gene_id ",'"',g,'"; transcript_id "',ex.names[t],'"; gene_name ','"',g,'"; gene_source "MEB"; gene_biotype "wes";',sep=""))
      x=paste(x,collapse = "\t")
      return(x)
    })    
    names(exGtf)=NULL
    write(exGtf, file=gtfoutFn, append = TRUE, sep="\t")
    #export exon
    exGtf=sapply(c(1:length(startR)),function(t){
      x=c(chr, "processed_transcript","exon",startR[t], endR[t],".",GEmat$EXONSTRAND[1],".", paste("gene_id ",'"',g,'"; transcript_id "',ex.names[t],'"; gene_name ','"',g,'"; gene_source "MEB"; gene_biotype "wes";',sep=""))
      x=paste(x,collapse = "\t")
      return(x)
    })    
    names(exGtf)=NULL
    write(exGtf, file=gtfoutFn, append = TRUE, sep="\t")
  }
}


rm(list=ls())

