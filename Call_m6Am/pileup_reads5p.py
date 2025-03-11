#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
import polars as pl
import argparse, time, pysam
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
    Pileup 5' ends of reads:

    1. reads with soft clip in 5' were excluded
    2. reads with too many unconverted As (>=3, not including the first A in read 5' end) were excluded
    3. use all A in downstream <window_size> bp or next A (if no A found in this window) as control
    """

    parser = argparse.ArgumentParser(prog="pileup_reads_5p",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-r","--ref", dest="references", nargs="+", required=True,help="reference fasta(s)")
    group_required.add_argument("-l","--list", dest="fTSS",required=True,
                                help="\nTSS list file: *_TSS_raw.bed.annotated.rmdup from anno_TSS.py. \
                                      \nFormat: ID,Chr,Start,End,Counts,Strand,Base,TPM,txChr,txStart,txEnd,geneID,txID,txStrand,txTSS,geneBiotype,txBiotype,Priorities,Dist,absDist,...; separated by tab.")
    group_required.add_argument("-b","--bam", dest="fbams", nargs="+", required=True,help="input bam(s), sorted")
    group_required.add_argument("-o","--output", dest="output",required=True,help="output")
    # Optional
    group_optional = parser.add_argument_group("Optional")
    group_optional.add_argument("-w", "--window", dest="window_size", default=30,
                                help="downstream window size to count unconverted As (as control), default=30. If no A found in this window, then use the first A found downstream as control.")
    group_optional.add_argument("--maxAs", dest="max_allowed_As", default=3,
                            help="maximum allowed number of As in reads to be count as signal ones, default=3")
    options = parser.parse_args()

    # Step1: Init
    ## parse reference
    print("------ [%s] Loading genome ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    reference_genome = {}
    for reference in options.references:
        for seq in SeqIO.parse(reference,"fasta"):
            reference_genome[seq.id] = str(seq.seq).upper()

    ## init site dict
    print("------ [%s] Initializing dict ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    output = {}
    with open(options.fTSS, "r") as fTSS:
        for line in fTSS:
            line = line.strip().split("\t")
            chr, pos, strand, base = line[1], line[3], line[5], line[6]
            if base.upper() != "A":
                continue

            pos = int(pos)
            ID = (chr, pos, strand)

            output[ID] = {"Chr": chr, "Pos": pos, "Strand": strand, "Base": base,  
                          "geneID": line[11], "txID": line[12], "txBiotype": line[16], "Dist": line[18],
                          "Counts": line[4], "TPM": line[7], "A": 0, "T": 0, "C": 0, "G": 0}
    
    ## find the position of next A 
    print("------ [%s] Finding next As ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    with open(options.fTSS, "r") as fTSS:
        for line in fTSS:
            line = line.strip().split("\t")
            chr, pos, strand, base = line[1], line[3], line[5], line[6]
            if base.upper() != "A":
                continue

            pos = int(pos)
            ID = (chr, pos, strand)

            # find next A
            lnext = "" # a list of next pos
            next_pos = pos
            idx = 0
            glen = len(reference_genome[chr])

            if strand == "+":
                while idx < options.window_size or not lnext:
                    next_pos += 1
                    if next_pos >= glen:
                        break
                    idx += 1
                    #key = (chr, next_pos, strand)
                    # if key not in output:
                    next_base = reference_genome[chr][next_pos-1]
                    if next_base == "A":
                        lnext = lnext + str(next_pos) + ","
            else:
                while idx < options.window_size or not lnext:
                    next_pos -= 1
                    if next_pos < 0:
                        break
                    idx += 1
                    #key = (chr, next_pos, strand)
                    # if key not in output:
                    next_base = reference_genome[chr][next_pos-1]
                    if next_base == "T":
                        lnext = lnext + str(next_pos) + ","

            output[ID]["Next_pos_A"] = 0
            output[ID]["Next_pos_T"] = 0
            output[ID]["Next_pos_C"] = 0
            output[ID]["Next_pos_G"] = 0
            output[ID]["Next_pos"] = lnext

    # Step II: Pileup 
    print("------ [%s] Pileup 5p ends of reads ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    for fbam in options.fbams:
        with pysam.AlignmentFile(fbam, "rb") as input_BAM:
            # iterate by reads, pileup 5' end
            for read in input_BAM:
                aligned_pairs = read.get_aligned_pairs()
                aligned_query_pos = [i[0] for i in aligned_pairs]
                aligned_ref_pos = [i[1] for i in aligned_pairs] # 0-based reference pos

                # get aligned positon of the 5' end of reads
                if read.is_forward:
                    pos = aligned_ref_pos[0]
                    query_base = read.query_sequence[0]
                    A_counts = read.query_sequence.count('A') # number of un-converted As in reads
                    if query_base == "A":
                        A_counts -= 1 # exclude the first A in reads
                else:
                    pos = aligned_ref_pos[-1]
                    query_base = reverse_complement(read.query_sequence[-1])
                    A_counts = read.query_sequence.count('T')
                    if query_base == "T":
                        A_counts -= 1

                if pos is not None: # exclude 5' soft clipped reads
                    pos = pos + 1 # convert to 1-based reference pos
                    chr = read.reference_name.replace("_AG_converted", "") # remove suffix
                    strand = '+' if read.is_forward else '-'
                    ID = (chr, pos, strand)

                    if ID in output and query_base in ("A", "T", "C", "G"):
                        if A_counts <= options.max_allowed_As: # only count signal reads (un-converted As <= max_allowed_As)
                            output[ID][query_base] += 1

                            # if next A pos is covered by this read, then count
                            lnext = output[ID]["Next_pos"].strip(",").split(",")
                            for next_pos in lnext:
                                try: 
                                    next_query_idx = aligned_ref_pos.index(int(next_pos)-1)
                                except ValueError:
                                    next_query_idx = None
                                
                                if next_query_idx is not None:
                                    next_query_pos = aligned_query_pos[next_query_idx]
                                    if next_query_pos is not None:
                                        next_query_base = read.query_sequence[next_query_pos]
                                        if read.is_reverse:
                                            next_query_base = reverse_complement(next_query_base)
                                        if next_query_base in ("A", "T", "C", "G"):
                                            output[ID]["Next_pos_"+next_query_base] += 1
                    
    # reformat data columns and save to file
    df = pl.from_pandas(pd.DataFrame.from_dict(output, orient='index'))
    info_col = ["Chr", "Pos", "Strand", "Base", "geneID", "txID", "txBiotype", "Dist"]
    data_col = ["A", "T", "C", "G", "Next_pos_A", "Next_pos_T", "Next_pos_C", "Next_pos_G", "Next_pos"]
    df.select(
        pl.col(info_col+data_col)
    ).write_csv(options.output, separator='\t')