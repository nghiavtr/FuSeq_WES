#!/usr/bin/env python

import os
import argparse
import re
import json
import pysam


cigarPattern = '([0-9]+[MIDNSHP])'
cigarSearch = re.compile(cigarPattern)
atomicCigarPattern = '([0-9]+)([MIDNSHP])'
atomicCigarSearch = re.compile(atomicCigarPattern)


def load_gtf(genes_gtf):
    # one time process, converting gtf file to json with transcripts, start, end
    gtf_lookup = dict()
    with open(genes_gtf, 'r') as fh:
        for line in fh.readlines():
            if '\ttranscript\t' in line:
                data = line.strip().split('\t')
                chrom = data[0].replace('chr', '')
                start = int(data[3])
                end = int(data[4])
                strand = data[6]
                info = data[8].split(';')
                transcript_id = info[1].split()[1].replace('"', '')
                gene_id = info[0].split()[1].replace('"', '')
                if '-AS1' in gene_id:
                    for key, pos in gtf_lookup[chrom].items():
                        gene_name = get_genename(key)
                        if gene_name in gene_id:
                            if (pos[0] <= start and start <= pos[1]) or \
                                (pos[0] <= end and end <= pos[1]):
                                gtf_lookup[chrom][key].append(transcript_id)
                else:
                    if chrom in gtf_lookup:
                        gtf_lookup[chrom].update({transcript_id: [start, end]})
                    else:
                        gtf_lookup[chrom] = {transcript_id: [start, end]}
    return gtf_lookup


def extractCigarOps(cigar, flag):
    if (cigar == "*"):
        cigarOps = []
    elif (flag & 0x0010):
        cigarOpStrings = cigarSearch.findall(cigar)
        cigarOps = []
        for opString in cigarOpStrings:
            cigarOpList = atomicCigarSearch.findall(opString)
            # "struct" for the op and it's length
            cigar = cigarOp(cigarOpList[0][0], cigarOpList[0][1])
            # add to the list of cigarOps
            cigarOps.append(cigar)
            cigarOps = cigarOps
        cigarOps.reverse()
        ##do in reverse order because negative strand##
    else:
        cigarOpStrings = cigarSearch.findall(cigar)
        cigarOps = []
        for opString in cigarOpStrings:
            cigarOpList = atomicCigarSearch.findall(opString)
            # "struct" for the op and it's length
            cigar = cigarOp(cigarOpList[0][0], cigarOpList[0][1])
            # add to the list of cigarOps
            cigarOps.append(cigar)
            # cigarOps = cigarOps
    return(cigarOps)


def calcQueryPosFromCigar(cigarOps):
    qsPos = 0
    qePos = 0
    qLen = 0
    # if first op is a H, need to shift start position
    # the opPosition counter sees if the for loop is looking at the first index of the cigar object
    opPosition = 0
    for cigar in cigarOps:
        if opPosition == 0 and (cigar.op == 'H' or cigar.op == 'S'):
            qsPos += cigar.length
            qePos += cigar.length
            #qLen += cigar.length
        elif opPosition > 0 and (cigar.op == 'H' or cigar.op == 'S'):
            qLen += 0
        elif cigar.op == 'M' or cigar.op == 'I':
            qePos += cigar.length
            qLen += cigar.length
            opPosition += 1
    d = queryPos(qsPos, qePos, qLen)
    return d


class cigarOp (object):
    """
    sturct to store a discrete CIGAR operations
    """
    def __init__(self, opLength, op):
        self.length = int(opLength)
        self.op     = op


class queryPos (object):
    """
    struct to store the start and end positions of query CIGAR operations
    """
    def __init__(self, qsPos, qePos, qLen):
        self.qsPos = int(qsPos)
        self.qePos = int(qePos)
        self.qLen  = int(qLen)


def calcQueryOverlap(s1, e1, s2, e2):
    o = 1 + min(e1, e2) - max(s1, s2)
    return max(0, o)


def checkOverlap(readCigar, mateCigar, readFlag, mateFlag, minNonOverlap):
    readCigarOps = extractCigarOps(readCigar, readFlag)
    readQueryPos = calcQueryPosFromCigar(readCigarOps)
    mateCigarOps = extractCigarOps(mateCigar, mateFlag)
    mateQueryPos = calcQueryPosFromCigar(mateCigarOps)
    overlap = calcQueryOverlap(readQueryPos.qsPos, readQueryPos.qePos, mateQueryPos.qsPos, mateQueryPos.qePos)
    nonOverlap1 = 1 + readQueryPos.qePos - readQueryPos.qsPos - overlap
    nonOverlap2 = 1 + mateQueryPos.qePos - mateQueryPos.qsPos - overlap
    mno = min(nonOverlap1, nonOverlap2)
    if mno >= minNonOverlap:
        return True
    return False


def get_transcript(chrom, start):
    if chrom in gtf_lookup:
        for transcript, pos in gtf_lookup[chrom].items():
            if pos[0] <= start and start <= int(pos[1]):
                if len(pos) > 2:
                    return [transcript, pos[2]]
                else:
                    return [transcript]
    return ''


def get_genename(transcript):
    if transcript:
        return transcript.split('__')[0]
    return ''


def calcHitpos(readpos, txts):
    if len(txts) == 1:
        return [int(readpos) - int(txts[0].split('__')[2])]
    else:
        hitpos = []
        for txt in txts:
            hitpos.append(int(readpos) - int(txt.split('__')[2]))
        return hitpos


def check_feq(discordant_read, feq_class):
    for ds_read, count in feq_class.items():
        if ds_read.ref_txts == discordant_read.ref_txts and ds_read.mate_txts == discordant_read.mate_txts:
            return ds_read, True
    return discordant_read, False


class discordantReads(object):
    def __init__(self, read, frst_txts, mate_txts):
        self.ref_txts = ','.join(frst_txts)
        self.mate_txts = ','.join(mate_txts)
        self.ref_rtype = 1 if read.is_read1 else 2
        self.mate_rtype = 2 if read.is_read1 else 1


parser = argparse.ArgumentParser(description="Extract discordant, split reads and build fuseq classes")
parser.add_argument("--bam", required=True, help="Input bam")
parser.add_argument("--gtf", required=True, help="Genes annotation (.json) file")
parser.add_argument("--mapq-filter", action='store_true', help="filter reads with Mapq<30")
parser.add_argument("--numSplits", default=2, help="The maximum number of split-read \
                    mappings to allow per read. Reads with more are excluded. Default=2", type=int)
parser.add_argument("--minNonOverlap", default=20, help="minimum non-overlap\
                    between split alignments on the query (default=20)", type=int)
parser.add_argument("--outdir", required=True, help="output dir")

args = parser.parse_args()
numSplits = args.numSplits
minNonOverlap = args.minNonOverlap
outdir = args.outdir
mapq = 30 if args.mapq_filter else 0

sr_output = open(os.path.join(outdir, "splitReadInfo.txt"), 'w')
mr_feq = open(os.path.join(outdir, "feq_ALL.txt"), 'w')

samfile = pysam.AlignmentFile(args.bam, "rb", threads=8)
reads = set()
read_ftx = dict()
feq_class = dict()
split_reads = set()

with open(args.gtf, 'r') as jsonfile:
    gtf_lookup = json.loads(jsonfile.read())

# gtf_json = open("/proj/snic2020-6-4/sarath/Results/BeatAML/UCSC_hg19_wes_contigSize3000_bigLen130000_r76.json", "w")
# gtf_json.write(json.dumps(gtf_lookup, indent = 4))

# processing bam file line by line
for read in samfile.fetch():
    # 1. Split Extraction start
    # split read pairs reported in three lines in bam file
    # but we can extract all information from one entry 
    split = False

    read_chrom = samfile.get_reference_name(read.reference_id).replace('chr', '')
    mate_chrom = read_chrom

    if samfile.get_reference_name(read.next_reference_id):
        mate_chrom = samfile.get_reference_name(read.next_reference_id).replace('chr', '')

    # print(read_chrom, mate_chrom)
    # break

    # to avoid repeat in split read entry
    if read.query_name in split_reads:
        continue

    # extract split reads
    if read.mapping_quality >= mapq:
        for tag in read.tags:
            if "SA" in tag:
                if(len(tag[1].split(";"))) <= numSplits:
                    sa = tag[1].split(",")
                    saCigar = sa[3]
                    saFlag = int(0)
                    saMapq = int(sa[4])
                    if sa[2] == "-":
                        saFlag = int(16)
                    if saMapq >= mapq:
                        readCigarOps = extractCigarOps(read.cigarstring, read.flag)
                        readQueryPos = calcQueryPosFromCigar(readCigarOps)
                        mateCigarOps = extractCigarOps(saCigar, saFlag)
                        mateQueryPos = calcQueryPosFromCigar(mateCigarOps)
                        overlap = calcQueryOverlap(readQueryPos.qsPos, readQueryPos.qePos, mateQueryPos.qsPos, mateQueryPos.qePos)
                        nonOverlap1 = 1 + readQueryPos.qePos - readQueryPos.qsPos - overlap
                        nonOverlap2 = 1 + mateQueryPos.qePos - mateQueryPos.qsPos - overlap
                        mno = min(nonOverlap1, nonOverlap2)
                        if mno >= minNonOverlap:
                            split = True
                        # split = checkOverlap(read.cigarstring, saCigar, read.flag, saFlag, minNonOverlap)

    # Annotating split reads with transcript information
    if split:

        frst_txts = get_transcript(read_chrom, read.reference_start)
        mate_txts = get_transcript(mate_chrom, read.next_reference_start)

        saChrom = sa[0].replace('chr', '')
        saStart = sa[1]
        saStrand = sa[2]

        sec_txts = get_transcript(saChrom, int(saStart))

        frst_gene = ','.join([get_genename(frst_txt) for frst_txt in frst_txts])
        sec_gene = ','.join([get_genename(sec_txt) for sec_txt in sec_txts])
        mate_gene = ','.join([get_genename(mate_txt) for mate_txt in mate_txts])

        read_type = 1 if read.is_read1 else 2 
        read_direction = 1 if read.is_reverse else 0
        mate_direction = 0 if read.is_reverse else 1

        frst_hitpos = ','.join(map(str, calcHitpos(read.pos, frst_txts)))
        sec_hitpos = ','.join(map(str, calcHitpos(saStart, sec_txts)))
        mate_hitpos = ','.join(map(str, calcHitpos(read.mpos, mate_txts)))

        # Writing split reads into output file
        if frst_txts != '' and sec_txts != '' and mate_txts != '' and\
           frst_gene != sec_gene:
            sr_output.write("\t".join(map(str, [read.query_name, read_type, read_direction, ','.join(frst_txts), frst_gene,
                            frst_hitpos, readQueryPos.qsPos, readQueryPos.qLen, ','.join(sec_txts), sec_gene, sec_hitpos,
                            mateQueryPos.qsPos, mateQueryPos.qLen, mate_gene, mate_direction, mate_hitpos])) + '\n')

        split_reads.add(read.query_name)

    # 2. Discordant reads extraction
    # In bam file, discordant reads are denoted in two lines, ref read and mate read
    if read.mapping_quality >= mapq:
        if read.query_name in read_ftx:
            feq_dsread, feq_check = check_feq(read_ftx[read.query_name], feq_class)
            if feq_check:
                feq_class[feq_dsread] += 1
            else:
                feq_class[read_ftx[read.query_name]] = 1
            continue

    # Using flag value we can extract discordant reads
    if not (1294 & read.flag) and not split:
        # discordant_reads.add(read.query_name)

        ref_txts = get_transcript(read_chrom, read.reference_start)
        mate_txts = get_transcript(mate_chrom, read.next_reference_start)

        if ref_txts != '' and mate_txts != '' and ref_txts != mate_txts and\
           read.mapping_quality >= 30:
            read_ftx[read.query_name] = discordantReads(read, ref_txts, mate_txts)

            if read.query_name in reads:
                # if reads[read.query_name][0].mapping_quality >= 30:
                feq_dsread, feq_check = check_feq(read_ftx[read.query_name], feq_class)
                if feq_check:
                    feq_class[feq_dsread] += 1
                else:
                    feq_class[read_ftx[read.query_name]] = 1

    reads.add(read.query_name)

# Extracting fusion eq class
mr_feq.write('\t'.join(['Transcript', 'Count', 'Read', 'Feq\n']))
eq_classID = 1
for feq, count in feq_class.items():

    if ',' in feq.ref_txts:
        txts = feq.ref_txts.split(',')
        for txt in txts:
            mr_feq.write('\t'.join(map(str, [txt, count, feq.ref_rtype, eq_classID])) + '\n')
    else:
        mr_feq.write('\t'.join(map(str, [feq.ref_txts, count, feq.ref_rtype, eq_classID])) + '\n')

    if ',' in feq.mate_txts:
        multi_txt = True
        txts = feq.mate_txts.split(',')
        for txt in txts:
            mr_feq.write('\t'.join(map(str, [txt, count, feq.mate_rtype, eq_classID])) + '\n')
    else:
        mr_feq.write('\t'.join(map(str, [feq.mate_txts, count, feq.mate_rtype, eq_classID])) + '\n')

    eq_classID += 1
