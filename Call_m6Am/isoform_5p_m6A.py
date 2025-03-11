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
    Count all unconverted A of each A sites in 5' regions downstream of each TSSs.

    1. Use TSS list as input, 5' regions were defined as downstream <window_size> bp from TSSs
    2. Reads with soft clipped in 5' end and reads with too many unconverted As (>=3, not including the first A in read 5' end) were excluded
    3. Only reads starting from specific TSS were considered (eg. isform level)
    """

    parser = argparse.ArgumentParser(prog="isoform_5p_m6A",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-r","--ref", dest="references", nargs="+", required=True,help="reference fasta(s)")
    group_required.add_argument("-l","--list", dest="fTSS",required=True,
                                help="\nTSS list file: usually *.TSS.passed of each sample from merge_multi_sample.py. \
                                      \nFormat: Chr,Pos,Strand,Base,geneID,txID,txBiotype,Dist,Counts,TPM,...; separated by tab.")
    group_required.add_argument("-b","--bam", dest="fbams", nargs="+", required=True,help="input bam(s), sorted")
    group_required.add_argument("-o","--output", dest="output",required=True,help="output")
    # Optional
    group_optional = parser.add_argument_group("Optional")
    group_optional.add_argument("-w", "--window", type=int, dest="window_size", default=100,
                                help="downstream window size to count unconverted As (as control), default=100.")
    group_optional.add_argument("--maxAs", dest="max_allowed_As", type=int, default=3,
                            help="maximum allowed number of As in reads to be count as signal ones, default=3")
    options = parser.parse_args()

    class NextPos(dict):
        def __init__(self):
            super().__init__()  # Initialize as a dictionary
            self.update({"A": 0, "T": 0, "C": 0, "G": 0})  # Add default values

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
        next(fTSS)
        for line in fTSS:
            line = line.strip().split("\t")
            chr, pos, strand, base = line[0], line[1], line[2], line[3]

            pos = int(pos)
            ID = (chr, pos, strand)

            output[ID] = {"Chr": chr, "Pos": pos, "Strand": strand, "Base": base,
                          "geneID": line[4], "txID": line[5], "txBiotype": line[6], "Dist": line[7], 
                          "Counts": line[8], "TPM": line[9],
                          "A": 0, "T": 0, "C": 0, "G": 0}
    
    ## find downstream A sites
    print("------ [%s] Finding downstream A sites ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    with open(options.fTSS, "r") as fTSS:
        next(fTSS)
        for line in fTSS:
            line = line.strip().split("\t")
            chr, pos, strand, base = line[0], line[1], line[2], line[3]

            pos = int(pos)
            ID = (chr, pos, strand)

            # find next A
            dnext = {} # a dict of next pos
            next_pos = pos
            idx = 0
            glen = len(reference_genome[chr])

            if strand == "+":
                while idx < options.window_size:
                    next_pos += 1
                    if next_pos >= glen:
                        break
                    idx += 1
                    #key = (chr, next_pos, strand)
                    # if key not in output:
                    next_base = reference_genome[chr][next_pos-1]
                    if next_base == "A":
                        dnext[next_pos] = NextPos()
            else:
                while idx < options.window_size:
                    next_pos -= 1
                    if next_pos < 0:
                        break
                    idx += 1
                    #key = (chr, next_pos, strand)
                    # if key not in output:
                    next_base = reference_genome[chr][next_pos-1]
                    if next_base == "T":
                        dnext[next_pos] = NextPos()

            output[ID]["Next_pos"] = dnext

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
                            # next_pos is reference 1-based pos, next_query_pos is query 0-based pos
                            dnext = output[ID]["Next_pos"]
                            for next_pos in dnext:
                                try: 
                                    next_query_idx = aligned_ref_pos.index(next_pos-1)
                                except ValueError:
                                    next_query_idx = None
                                
                                if next_query_idx is not None:
                                    next_query_pos = aligned_query_pos[next_query_idx]
                                    if next_query_pos is not None:
                                        next_query_base = read.query_sequence[next_query_pos]
                                        if read.is_reverse:
                                            next_query_base = reverse_complement(next_query_base)
                                        if next_query_base in ("A", "T", "C", "G"):
                                            output[ID]["Next_pos"][next_pos][next_query_base] += 1
                    

with open(options.output, "w") as fout:
    tmp = next(iter(output))
    col_names = list(output[tmp].keys())
    fout.write("\t".join(col_names))
    fout.write("\t"+"Next_pos_ATCG"+"\n")

    for ID in output:
        info = [str(output[ID][col]) for col in col_names[:-1]]
        fout.write("\t".join(info)+"\t")  # information

        dnext = output[ID]["Next_pos"]
        fout.write(";".join([str(next_pos) for next_pos in dnext])+";\t")
        for next_pos in dnext:
            fout.write(",".join(
                [str(dnext[next_pos]['A']), str(dnext[next_pos]['T']), 
                 str(dnext[next_pos]['C']), str(dnext[next_pos]['G'])]) + ";")  # SP=bins[i]['total']/counts['total']
        fout.write("\n")

print("[%s] Done!" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)