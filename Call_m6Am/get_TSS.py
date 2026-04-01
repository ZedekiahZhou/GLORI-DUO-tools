#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
import argparse
import pysam
import pandas as pd
from Bio.Seq import reverse_complement

def get_refer_base(key):
    global reference_genome
    chr, pos, strand = key
    if strand == "+":
        return reference_genome[chr][pos - 1]
    else:
        return reverse_complement(reference_genome[chr][pos -1])


if __name__ == "__main__":
    description = """
    get transcription start sites by pile up 5' ends of non-softclipped reads.
    """

    parser = argparse.ArgumentParser(prog="get_TSS",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-r","--refs", dest="references", nargs='+', required=True, help="reference fasta(s)")
    group_required.add_argument("-b","--bam", dest="fbams", nargs='+', required=True, 
                                help="one or more sorted input bam (Multimapped and unmapped reads must be excluded first!)")
    group_required.add_argument("-o","--output",dest="output",required=True,help="output")
    options = parser.parse_args()

    # Step1: Init
    ## parse reference
    reference_genome = {}
    for reference in options.references:
        for seq in SeqIO.parse(reference,"fasta"):
            reference_genome[seq.id] = str(seq.seq).upper()

    ## init site dict
    output = {}
    total = 0
    mapped = 0

    # Pileup 
    for fbam in options.fbams:
        with pysam.AlignmentFile(fbam, "rb") as input_BAM:
            # iterate by reads, pileup 5' end
            for read in input_BAM:
                aligned_ref_pos = read.get_reference_positions(full_length=True) # 0-based reference pos
                pos = aligned_ref_pos[-1] if read.is_reverse else aligned_ref_pos[0] # get aligned positon of the 5' end of reads

                if pos is not None: # exclude 5' soft clipped reads
                    pos = pos + 1 # convert to 1-based reference pos
                    chr = read.reference_name.replace("_AG_converted", "") # remove suffix
                    strand = '-' if read.is_reverse else '+'
                    ID = (chr, pos, strand)
                    if ID not in output:
                        output[ID] = {}
                        output[ID]["Counts"] = 0
                        output[ID]["Ref_base"] = get_refer_base(ID)
                    
                    output[ID]["Counts"] += 1
                    total += 1
                mapped += 1
        
    # get TPM
    with open(options.output, "w") as fout:
        for key, values in output.items():
            ID = key[0] + "_" + str(key[1]) + "_" + key[2]
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key[0], key[1]-1, key[1], ID, values["Counts"], key[2], 
                                                                 values["Ref_base"], round(values["Counts"]/total*1000000, 3)))

    print("Non-softclipped reads: {}, Mapped reads: {}, Non-softclipped rate: {}%".
          format(total, mapped, round(total/mapped*100, 3)))    