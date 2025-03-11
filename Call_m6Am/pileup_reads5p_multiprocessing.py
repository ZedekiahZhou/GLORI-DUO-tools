#!/usr/bin/env python

from Bio import SeqIO
import argparse, time, pysam
import pandas as pd
import polars as pl
from Bio.Seq import reverse_complement
import multiprocessing


def get_refer_base(key):
    global reference_genome
    chr, pos, strand = key
    if strand == "+":
        return reference_genome[chr][pos - 1]
    else:
        return reverse_complement(reference_genome[chr][pos -1])


class NextPos(dict):
    def __init__(self):
        super().__init__()  # Initialize as a dictionary
        self.update({"A": 0, "T": 0, "C": 0, "G": 0})  # Add default values


def pileup_bin(chr, bin_start, bin_end, bin_df):
    global options
    print("[%s] Pileup bin %s:%d-%d" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), 
                                        bin_chr, bin_start, bin_end), flush=True)
    if bin_df.shape[0] == 0:
        return None
    else:
        res = {}
        for row in bin_df.iter_rows():
            ID = (row[0], row[1], row[2])
            res[ID] = {"Chr": row[0], "Pos": row[1], "Strand": row[2], "Base": row[3],  
                       "geneID": row[4], "txID": row[5], "txBiotype": row[6], "Dist": row[7],
                       "Counts": row[8], "TPM": row[9], 
                       "A": 0, "T": 0, "C": 0, "G": 0}
            
            res[ID]["Next_pos"] = {} # a dict of next pos
            for next_pos in row[10]:
                res[ID]["Next_pos"][next_pos] = NextPos()

    for fbam in options.fbams:
        try:
            with pysam.AlignmentFile(fbam, "rb") as input_BAM:
                # iterate by reads, pileup 5' end
                for read in input_BAM.fetch(contig=bin_chr+"_AG_converted", start=bin_start, end=bin_end):
                    read5p = read.reference_start if read.is_forward else read.reference_end-1 # 0-based reference pos of 5' end of reads
                    if not (read5p >= bin_start and read5p < bin_end):
                        continue # skip reads with 5' end outside of bin
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

                        if ID in res and query_base in ("A", "T", "C", "G"):
                            if A_counts <= options.max_allowed_As: # only count signal reads (un-converted As <= max_allowed_As)
                                res[ID][query_base] += 1

                                # if next A pos is covered by this read, then count
                                # next_pos is reference 1-based pos, next_query_pos is query 0-based pos
                                for next_pos in res[ID]["Next_pos"]:
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
                                                res[ID]["Next_pos"][next_pos][next_query_base] += 1
        except ValueError:
            print("[%s] No bin %s:%d-%d in %s" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), 
                                        bin_chr, bin_start, bin_end, fbam), flush=True)
    return res


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
    group_optional.add_argument("-s", "--step", dest="step_size", type=int, default=500000, 
                                help="step size to process bam file, default=500000.")
    group_optional.add_argument("-p", "--processes", dest="processes", type=int, default=4,  
                                help="number of processes to use, default=4.")
    group_optional.add_argument("--maxAs", dest="max_allowed_As", default=3,
                            help="maximum allowed number of As in reads to be count as signal ones, default=3")
    options = parser.parse_args()

    # Step1: Init
    ## parse reference
    print("------ [%s] Loading genome ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    reference_genome = {}
    RefBins = []
    for reference in options.references:
        for seq in SeqIO.parse(reference,"fasta"):
            reference_genome[seq.id] = str(seq.seq).upper()
            for bin_start in range(0, len(seq.seq), options.step_size):  # split reference into bins
                RefBins.append((seq.id, bin_start, bin_start + options.step_size))
    
    ## init site dict
    print("------ [%s] Initializing dict ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    tss = dict()
    with open(options.fTSS, "r") as fTSS:
        for line in fTSS:
            line = line.strip().split("\t")
            chr, pos, strand, base = line[1], line[3], line[5], line[6]
            if base.upper() != "A":
                continue

            pos = int(pos)
            ID = (chr, pos, strand)

            tss[ID] = {"Chr": chr, "Pos": pos, "Strand": strand, "Base": base,  
                          "geneID": line[11], "txID": line[12], "txBiotype": line[16], "Dist": line[18],
                          "Counts": line[4], "TPM": line[7]}
    
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
            dnext = []
            next_pos = pos
            idx = 0
            glen = len(reference_genome[chr])

            if strand == "+":
                while idx < options.window_size or not dnext:
                    next_pos += 1
                    if next_pos >= glen:
                        break
                    idx += 1
                    #key = (chr, next_pos, strand)
                    # if key not in output:
                    next_base = reference_genome[chr][next_pos-1]
                    if next_base == "A":
                        dnext.append(next_pos)
            else:
                while idx < options.window_size or not dnext:
                    next_pos -= 1
                    if next_pos < 0:
                        break
                    idx += 1
                    #key = (chr, next_pos, strand)
                    # if key not in output:
                    next_base = reference_genome[chr][next_pos-1]
                    if next_base == "T":
                        dnext.append(next_pos)

            tss[ID]["Next_pos"] = dnext

    df = pl.from_pandas(pd.DataFrame.from_dict(tss, orient='index'))

    # Step II: Pileup 
    print("------ [%s] Pileup 5p ends of reads ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    pool = multiprocessing.Pool(processes=options.processes)

    async_results = []
    t1=time.time()
    try:
        for item in RefBins:
            bin_chr, bin_start, bin_end = item
            bin_df = df.filter((pl.col("Chr") == bin_chr) & (pl.col("Pos") > bin_start) & 
                               (pl.col("Pos") <= bin_end)) # Pos is 1-based, bin is 0-based
            async_result = pool.apply_async(pileup_bin, args=(bin_chr, bin_start, bin_end, bin_df))
            async_results.append(async_result)

        pool.close()
        print("------ [%s] Run multiprocessing using [%s] cores ..." % 
              (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), str(options.processes)), flush=True)
        pool.join()
        print("------ [%s] Join results ..." % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    finally:
        pool.terminate()

    t2=time.time()
    print("------ [%s] Pileup time used: %.2f s" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), t2-t1), flush=True)


    # write out
    with open(options.output, "w") as fout:
        col_names = ["Chr", "Pos", "Strand", "Base", "geneID", "txID", "txBiotype", "Dist", 
                     "Counts", "TPM", "A", "T", "C", "G"]
        fout.write("\t".join(col_names +["Next_pos_A", "Next_pos_T", "Next_pos_C", "Next_pos_G", "Next_pos"]) + "\n")
        for async_result in async_results:
            output = async_result.get()
            if output:
                for ID in output:
                    info = [str(output[ID][col]) for col in col_names]
                    fout.write("\t".join(info)+"\t")  # information

                    dnext = output[ID]["Next_pos"]
                    next_ATCG = [sum([dnext[next_pos][used_base] for next_pos in dnext]) for used_base in ["A", "T", "C", "G"]]
                    fout.write("\t".join(next_ATCG)+"\t") # sum of all next pos
                    fout.write(",".join([str(next_pos) for next_pos in dnext])+",\n")

    print("[%s] Done!" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)