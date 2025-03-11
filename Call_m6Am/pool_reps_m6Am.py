#!/usr/bin/env python
"""
Author: Zhe Zhou, Peking University, Yi lab
Date: Dec 27, 2024
Email: zzhou24@pku.edu.cn
Program: Merge replicates into one sample, recalculate TPM
used polars version: 0.20.25 --> 1.5.0
"""

import argparse, re, time
import polars as pl

parser = argparse.ArgumentParser(description="Pool TSS (and m6Am) sites from replicates")
parser.add_argument("-i", "--inputs", type=str, nargs='+', required=True, 
                    help="Input files recording TSS sites (*_TSS_raw.bed.annotated.rmdup)")
parser.add_argument("-o", "--output", type=str, required=True, 
                    help="output file prefix (output file will beprefix_TSS_raw.bed.annotated.rmdup and prefix_m6Am_sites_raw.tsv)")
parser.add_argument("--untreated", action="store_true", 
                    help="Untreated samaples (only merge TSS)")
args = parser.parse_args()


# ID      Chr     Start   End     Counts  Strand  Base    TPM     
# txChr   txStart txEnd   geneID  txID    txStrand        txTSS   
# geneBiotype     txBiotype       Priorities      Dist    absDist 
# raw_totalTPM    raw_nTSS        meanTPM stdTPM    relSum  zscore
tss_col = ["Chr", "End", "Strand", "Base",  
           "geneID", "txID", "txBiotype", "Dist", "absDist"]

# I. merge TSS sites, recalculate TPM
# read and merge data --------
print("Pool TSS sites: ", flush=True)
for i in range(len(args.inputs)):
    print("\t[%s] Read sample %s ========" % 
        (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), args.inputs[i]), flush=True)

    tss = pl.read_csv(args.inputs[i], separator="\t")

    # outer join samples
    if i == 0:
        tss_merged = tss
    else:
        tss_merged = tss_merged.join(
            tss, suffix="_right",
            on=tss_col, 
            how="full", coalesce=True
        ).with_columns(
            Counts=pl.sum_horizontal("Counts", "Counts_right"),
        ).select(tss_col + ["Counts"])

# calculate TPM, group sum, mean, and std
tss_merged = tss_merged.with_columns(
    TPM=pl.col("Counts") / pl.col("Counts").sum() * 1000000,
)

tss_agg=tss_merged.group_by("geneID").agg(
    pl.col("TPM").sum().alias("raw_totalTPM"),
    pl.col("TPM").count().alias("raw_nTSS"),
    pl.col("TPM").mean().alias("meanTPM"),
    pl.col("TPM").std().alias("stdTPM")
)

tss_merged = tss_merged.join(tss_agg, on="geneID", how="left").with_columns(
    (pl.col("TPM")/pl.col("raw_totalTPM")).round(3).alias("relSum"),
    pl.when((pl.col("raw_nTSS") <= 3))
    .then(pl.lit(99))
    .when(pl.col("stdTPM") == 0)
    .then(pl.lit(0))
    .otherwise(((pl.col("TPM") - pl.col("meanTPM")) / pl.col("stdTPM")))
    .round(3).alias("zscore"), 
    pl.col("raw_totalTPM").round(3),
    pl.col("meanTPM").round(3),
    pl.col("stdTPM").round(3)
)

tss_merged = tss_merged.sort(["Chr", "End"])
tss_merged.write_csv(args.output+"_TSS_raw.bed.annotated.rmdup", separator="\t",null_value=".")


# II. merged m6Am sites
# Chr     Pos     Strand  Base        Counts  TPM     geneID  txID    txBiotype       Dist
# A       T       C       G       Next_pos_A      Next_pos_T      Next_pos_C      Next_pos_G      Next_pos
# Signal_cov      AG_cov  Next_pos_AG     Signal_Ratio    AG_Ratio        Ctrl_Ratio      m6Am_Ratio      Pvalue  FDR
m6Am_cols = ["Chr", "Pos", "Strand", "Base", "geneID", "txID", "txBiotype", "Dist", 
             "AG_cov", "A", "Signal_Ratio", "AG_Ratio", "m6Am_Ratio", "Pvalue", "FDR"]
if not args.untreated:
    print("Pool m6Am sites: ", flush=True)
    for i in range(len(args.inputs)):
        tmp_input = args.inputs[i].replace("_TSS_raw.bed.annotated.rmdup", "_m6Am_sites_raw.tsv")

        print("\t[%s] Read sample %s ========" % 
            (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), tmp_input), flush=True)

        m6Am = pl.read_csv(tmp_input,separator="\t", null_values=".")

        # outer join samples, m6Am_Ratio not updated
        if i == 0:
            m6Am_merged = m6Am
        else:
            m6Am_merged = m6Am_merged.join(
                m6Am, suffix="_right",
                on=["Chr", "Pos", "Strand", "Base", "geneID", "txID", "txBiotype", "Dist"], 
                how="full", coalesce=True 
            ).with_columns(
                AG_cov=pl.sum_horizontal("AG_cov", "AG_cov_right"),
                A=pl.sum_horizontal("A", "A_right"),
                Signal_Ratio=pl.mean_horizontal("Signal_Ratio", "Signal_Ratio_right"),
                AG_Ratio=pl.mean_horizontal("AG_Ratio", "AG_Ratio_right"),
                Pvalue=pl.min_horizontal("Pvalue", "Pvalue_right"),
                FDR=pl.min_horizontal("FDR", "FDR_right"),
            ).select(m6Am_cols)
    
    m6Am_merged = m6Am_merged.with_columns(
        m6Am_Ratio=(pl.col("A") / pl.col("AG_cov")).round(3)
    ).select(m6Am_cols)
    m6Am_merged = m6Am_merged.sort(["Chr", "Pos"])
    m6Am_merged.write_csv(args.output+"_m6Am_sites_raw.tsv", separator="\t",null_value=".")