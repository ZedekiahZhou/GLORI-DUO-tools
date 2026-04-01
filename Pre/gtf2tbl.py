#!/usr/bin/env python

# Author: Zhou Zhe
# Date: Mar 25, 2024
# Email: cicadapoplar@163.com
# Usage: Extract GTF annotation to annotation table (UCSC table like)
# Input: [.gtf]

import argparse

if __name__ == "__main__":
	#Parser
    usage = "Usage: python %prog -i <gtf> -o output.tbl"
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input", required=True, help="Input GTF file, both NCBI and gencode GTFs were supported!")
    parser.add_argument("-o", dest="output", required=False, help="Output *.anno file")
    options = parser.parse_args()

    #format:
    # transcript_id chr strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score gene_id gene_name gene_biotype transcript_biotype transcript_length
    input = options.input
    output = options.output

    if not output:
        output = input[:-3] + "tbl"
    output_longest = output + ".longest"

    with open(input, "r") as gtf:
        line = gtf.readline()
        dict_gene = {}
        dict_trans = {}
        
        while line:
            if line.startswith("#"): 
                line = gtf.readline()
                continue
            
            # get basic info of lines
            row = line.strip().split(sep = "\t")
            row2 = row[8].split()
            rtype = row[2]
            gene_id = row2[row2.index("gene_id")+1].strip('";')
            try:
                transcript_id = row2[row2.index("transcript_id")+1].strip('";')
            except ValueError:
                transcript_id = "NA"
            
            if rtype == "gene":
                if gene_id not in dict_gene:
                    try:
                        gene_biotype = row2[row2.index("gene_biotype")+1].strip('";')
                    except ValueError:
                        try:
                            gene_biotype = row2[row2.index("gene_type")+1].strip('";')
                        except ValueError:
                            gene_biotype = "NA"
                    dict_gene[gene_id] = {"gene_biotype": gene_biotype, "transcripts": []}
                else:
                    print("Duplicate gene ID: ", gene_id)
            elif rtype == "transcript":
                if transcript_id not in dict_trans:
                    try:
                        transcript_biotype = row2[row2.index("transcript_biotype")+1].strip('";')
                    except ValueError:
                        try:
                            transcript_biotype = row2[row2.index("transcript_type")+1].strip('";')
                        except ValueError:
                            transcript_biotype = "NA"

                    try:
                        gene_name = row2[row2.index("gene_name")+1].strip('";')
                    except ValueError:
                        gene_name = gene_id

                    dict_trans[transcript_id] = {'trans_id': transcript_id, 'chr': row[0], 'strand': row[6], 'gene_id': gene_id,
                                                'txStart': int(row[3])-1, 'txEnd': int(row[4]), 
                                                'cdsStart': '.', 'cdsEnd': '.', 
                                                'exonStarts': [], 'exonEnds': [],
                                                'gene_id': gene_id, 'gene_name': gene_name,
                                                'transcript_biotype': transcript_biotype, 
                                                'gene_biotype': dict_gene[gene_id]["gene_biotype"]}
                    dict_gene[gene_id]["transcripts"] += [transcript_id]
                else:
                    print("Duplicate transcript ID: ", transcript_id)
            elif transcript_id not in dict_trans:  # unknown transcripts
                line = gtf.readline()
                continue
            elif rtype == "start_codon":  # cds 
                if row[6] == "+":
                    dict_trans[transcript_id]['cdsStart'] = int(row[3])-1
                else:
                    dict_trans[transcript_id]['cdsEnd'] = int(row[4])
            elif rtype == "stop_codon":  # cds 
                if row[6] == "+":
                    dict_trans[transcript_id]['cdsEnd'] = int(row[3])-1
                else:
                    dict_trans[transcript_id]['cdsStart'] = int(row[4])
            elif rtype == "exon":
                dict_trans[transcript_id]['exonStarts'] += [int(row[3])-1]
                dict_trans[transcript_id]['exonEnds'] += [int(row[4])]
            
            # get new line
            line = gtf.readline()

    # write out
    header = 'transcript_id chr strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score gene_id gene_name gene_biotype transcript_biotype transcript_length'.split()
    with open(output, "w") as anno:
        anno.writelines('\t'.join(header) + '\n')
        for trans_id in dict_trans.keys():
            trans = dict_trans[trans_id]
            trans['exonCount'] = len(trans['exonStarts'])
            exonStarts = ','.join([str(i) for i in sorted(trans['exonStarts'])]) + ","
            exonEnds = ','.join([str(i) for i in sorted(trans['exonEnds'])]) + ","
            trans['transcript_length'] = sum([trans['exonEnds'][i] - trans['exonStarts'][i] for i in range(trans['exonCount'])])
            anno.writelines('\t'.join([trans['trans_id'], trans['chr'], trans['strand'], 
                                       str(trans['txStart']), str(trans['txEnd']),
                                       str(trans['cdsStart']), str(trans['cdsEnd']), str(trans['exonCount']), exonStarts, exonEnds,
                                       ".", trans['gene_id'], trans['gene_name'],
                                       trans['gene_biotype'], trans['transcript_biotype'], str(trans['transcript_length'])]) + "\n")
    
    # get longest transcripts
    dict_longest = {}
    for trans_id in dict_trans.keys():
        trans = dict_trans[trans_id]
        if (trans['gene_id'] not in dict_longest) or (trans['transcript_length'] > dict_longest[trans['gene_id']]['transcript_length']):
            dict_longest[trans['gene_id']] = trans

    # write out
    with open(output_longest, "w") as anno:
        anno.writelines('\t'.join(header) + '\n')
        for gene_id in dict_longest.keys():
            trans = dict_longest[gene_id]
            exonStarts = ','.join([str(i) for i in sorted(trans['exonStarts'])]) + ","
            exonEnds = ','.join([str(i) for i in sorted(trans['exonEnds'])]) + ","
            anno.writelines('\t'.join([trans['trans_id'], trans['chr'], trans['strand'], 
                                       str(trans['txStart']), str(trans['txEnd']),
                                       str(trans['cdsStart']), str(trans['cdsEnd']), str(trans['exonCount']), exonStarts, exonEnds,
                                       ".", trans['gene_id'], trans['gene_name'],
                                       trans['gene_biotype'], trans['transcript_biotype'], str(trans['transcript_length'])]) + "\n")
    
    

