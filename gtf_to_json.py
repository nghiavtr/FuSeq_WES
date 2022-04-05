#!/usr/bin/env python

import argparse
import json


def load_gtf(genes_gtf):
    # one time process, converting gtf file to json with transcripts, start, end
    gtf_lookup = dict()
    with open(genes_gtf, 'r') as fh:
        for line in fh.readlines():
            if '\ttranscript\t' in line:
                data = line.strip().split('\t')
                chrom = data[0].replace('chr', '')
                if data[3] != 'NA':
                    start = int(data[3])
                    end = int(data[4])
                    strand = data[6]
                    info = data[8].split(';')
                    transcript_id = info[1].split()[1].replace('"', '')
                    gene_id = info[0].split()[1].replace('"', '')
					# skipping anti-strand genes, need to be fixed	
                    if '-AS1' in gene_id:
                        continue
                        # for key, pos in gtf_lookup[chrom].items():
                        #     gene_name = get_genename(key)
                        #     if gene_name in gene_id:
                        #         if (pos[0] <= start and start <= pos[1]) or \
                        #             (pos[0] <= end and end <= pos[1]):
                        #             gtf_lookup[chrom][key].append(transcript_id)
                    else:
                        if chrom in gtf_lookup:
                            gtf_lookup[chrom].update({transcript_id: [start, end]})
                        else:
                            gtf_lookup[chrom] = {transcript_id: [start, end]}
    return gtf_lookup


def get_genename(transcript):
    if transcript:
        return transcript.split('__')[0]
    return ''


parser = argparse.ArgumentParser(description="Convert gtf files to Json")
parser.add_argument("--gtf", required=True, help="Genes annotation file")
parser.add_argument("--output", required=True, help="Output json")

args = parser.parse_args()
gtf_file = args.gtf
outfname = args.output

gtf_lookup = load_gtf(gtf_file)

gtf_json = open(outfname, 'w')
gtf_json.write(json.dumps(gtf_lookup, indent = 4))
