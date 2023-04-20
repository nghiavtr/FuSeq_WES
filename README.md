# FuSeq_WES: a Fusion Detection method from DNA sequencing data

This tool is developed based on [FuSeq](https://github.com/nghiavtr/FuSeq), the method for detecting fusion genes from RNA-seq data. Many functions and parameters of FuSeq_WES are inherited from FuSeq. More details about FuSeq_WES can be found in [its article](https://www.frontiersin.org/article/10.3389/fgene.2022.820493).

## Installation Requirements

	- Git
	- Python3
	- pysam 
	- R (v3.6 or later) 
	- R packages - GenomicFeatures, biomaRt, Biostrings

## How to Install

Install pysam using pip

```sh
sudo apt update
sudo apt install python3-pip

pip install pysam
```
or 

Install pysam using conda

To install conda, please follow the instruction in this link https://docs.conda.io/projects/conda/en/latest/user-guide/install/ 

```sh
conda install -c bioconda pysam
```
## Quick start to run FuSeq-WES

```sh
# download FuSeq_WES
# If use FuSeq_WES_v1.0.0
wget https://github.com/nghiavtr/FuSeq_WES/releases/download/v1.0.0/FuSeq_WES_v1.0.0.tar.gz -O FuSeq_WES_v1.0.0.tar.gz
tar -xzvf FuSeq_WES_v1.0.0.tar.gz
mv FuSeq_WES_v1.0.0 FuSeq_WES

## If use the developing version from GitHub
#git clone https://github.com/nghiavtr/FuSeq_WES

#configure FuSeq_WES
cd FuSeq_WES
bash configure.sh
cd ..

# download test data
wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/2022/04/FuSeq_WES_testdata.tar.gz
tar -xzvf FuSeq_WES_testdata.tar.gz
rm FuSeq_WES_testdata.tar.gz

# download references
wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/2022/04/UCSC_hg19_wes_contigSize3000_bigLen130000_r100.tar.gz
tar -xzvf UCSC_hg19_wes_contigSize3000_bigLen130000_r100.tar.gz
rm UCSC_hg19_wes_contigSize3000_bigLen130000_r100.tar.gz
wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/2022/04/UCSC_hg38_wes_contigSize3000_bigLen130000_r100.tar.gz
tar -xzvf UCSC_hg38_wes_contigSize3000_bigLen130000_r100.tar.gz
rm UCSC_hg38_wes_contigSize3000_bigLen130000_r100.tar.gz

bamfile="FuSeq_WES_testdata/test.bam"
ref_json="UCSC_hg19_wes_contigSize3000_bigLen130000_r100/UCSC_hg19_wes_contigSize3000_bigLen130000_r100.json"
gtfSqlite="UCSC_hg19_wes_contigSize3000_bigLen130000_r100/UCSC_hg19_wes_contigSize3000_bigLen130000_r100.sqlite"

output_dir="test_out"
mkdir $output_dir

#extract mapped reads and split reads
python3 FuSeq_WES/fuseq_wes.py --bam $bamfile  --gtf $ref_json --mapq-filter --outdir $output_dir

#process the reads
fusiondbFn="FuSeq_WES/Data/Mitelman_fusiondb.RData"
paralogdbFn="FuSeq_WES/Data/ensmbl_paralogs_grch37.RData"
Rscript FuSeq_WES/process_fuseq_wes.R in=$output_dir sqlite=$gtfSqlite fusiondb=$fusiondbFn paralogdb=$paralogdbFn out=$output_dir

# Fusion genes discovered by FuSeq_WES are stored in a file named FuSeq_WES_FusionFinal.txt
# the other information of split reads and mapped reads are also founded in the output folder
```

## Building references for FuSeq_WES
This section shows an example of generating reference files (the json file and the sqlite file) for running fuseq_wes.
Users should select the right annotation version (hg19/hg38) and parameters (for example, the read length of bam files, default=100).

We also provide several pre-built references for FuSeq_WES that can be downloaded here: [Hg19 references](https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/2022/04/UCSC_hg19_wes_contigSize3000_bigLen130000_r100.tar.gz) and [Hg38 references](https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/2022/04/UCSC_hg38_wes_contigSize3000_bigLen130000_r100.tar.gz).

```sh
# download FuSeq_WES
# If use FuSeq_WES_v1.0.0
wget https://github.com/nghiavtr/FuSeq_WES/releases/download/v1.0.0/FuSeq_WES_v1.0.0.tar.gz -O FuSeq_WES_v1.0.0.tar.gz
tar -xzvf FuSeq_WES_v1.0.0.tar.gz
mv FuSeq_WES_v1.0.0 FuSeq_WES

## If use the developing version from GitHub
#git clone https://github.com/nghiavtr/FuSeq_WES

#configure FuSeq_WES
cd FuSeq_WES
bash configure.sh
cd ..

#download hg38 annotation from ucsc
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

gunzip hg38.fa.gz
gunzip hg38.refGene.gtf.gz

fusiondbFn="FuSeq_WES/Data/Mitelman_fusiondb.RData"
transcriptGtfSqlite="hg38.refGene.sqlite"
genomeFastaFile="hg38.fa"

#create sqlite file for the transcriptome
Rscript FuSeq_WES/createSqlite.R hg38.refGene.gtf $transcriptGtfSqlite

#create gtf file for all genes
readLen=100 #suppose the bam files have read length of 100bp. The the results can be slightly different if using this reference for input data with different read length.
gtfoutFn="UCSC_hg38_wes_r100.gtf"
Rscript FuSeq_WES/extract_gtf.R genomefasta=$genomeFastaFile sqlite=$transcriptGtfSqlite fusiondb=$fusiondbFn readLen=$readLen out=$gtfoutFn

### now generate references for FuSeq_WES
json="UCSC_hg38_wes_r100.json"
gtfSqliteOut="UCSC_hg38_wes_r100.sqlite"

#create sqlite file for the reference of FuSeq_WES
Rscript FuSeq_WES/createSqlite.R $gtfoutFn $gtfSqliteOut

#create JSON file for the reference of FuSeq_WES
python3 FuSeq_WES/gtf_to_json.py --gtf $gtfoutFn --output $json
``` 

#### Building the database of paralogs for FuSeq_WES
Information of paralogs is added to the output of FuSeq_WES, however this information is not used for filtering. Here, we show an example of scripts to get the data of paralogs for Hg38. Prebuilt data for paralog database of hg38 can be found in folder Data: https://github.com/nghiavtr/FuSeq_WES/tree/main/Data

1. Download gtf file from ensembl annotation database, create a sqlite file
``` 
#run in linux command line
wget http://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz
gunzip Homo_sapiens.GRCh38.95.gtf.gz
Rscript FuSeq_WES/createSqlite.R Homo_sapiens.GRCh38.95.gtf Homo_sapiens.GRCh38.95.sqlite
``` 

2. Extract paralog information and save to file
``` 
library(GenomicFeatures)
gtfSqlite=("Homo_sapiens.GRCh38.95.sqlite")
anntxdb <- loadDb(gtfSqlite)
allGenes=genes(anntxdb)
allGenes_name=names(allGenes)

library(biomaRt)
#remove version in the gene name before querying
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version = 95) #GRCh=38

att_meta=attributes(ensembl)
att_all=att_meta$attributes
head(att_all)

attList=c('ensembl_gene_id',"external_gene_name","hsapiens_paralog_ensembl_gene","hsapiens_paralog_associated_gene_name","hsapiens_paralog_orthology_type")
paralogs <- getBM(attributes=attList, filters = 'ensembl_gene_id', values = allGenes_name, mart = ensembl)

#rename the columnn
p=which(att_all$name %in% attList & att_all$page=="homologs")
newName=att_all$description[p]
newName=gsub(" ",".",newName)
colnames(paralogs)=newName

save(paralogs, file="ensmbl_paralogs_grch38.RData")
``` 

## References
1. Deng, Wenjiang, Sarath Murugan, Johan Lindberg, Venkatesh Chellappa, Xia Shen, Yudi Pawitan, and Trung Nghia Vu. 2022. “Fusion Gene Detection Using Whole-Exome Sequencing Data in Cancer Patients.” Frontiers in Genetics 13. https://www.frontiersin.org/article/10.3389/fgene.2022.820493.
