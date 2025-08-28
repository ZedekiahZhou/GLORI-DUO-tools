#!/usr/bin/env python
"""
Author: Zhe Zhou, Peking University, Yi lab
Date: Dec 27, 2024
Email: zzhou24@pku.edu.cn
Program: Merge and filter m6A sites from multiple samples
used polars version: 0.20.25 --> 1.5.0
"""

import argparse, re, time
import polars as pl
from pathlib import Path

parser = argparse.ArgumentParser(description="Merge and filter m6A sites from multiple samples")
parser.add_argument("-i", "--inputs", type=str, nargs='+', required=True, 
                    help="Input files recording m6A sites (*totalm6A.FDR.csv)")
parser.add_argument("--prx", type=str, nargs='+', 
                    help="columns names (usually sample names) to used for each input")
parser.add_argument("-o", "--output", type=str, required=True, help="output file name")
parser.add_argument("-c", "--AGcov", type=int, default=15, 
                    help='minimum A+G coverage for TSS (default: 15)')
parser.add_argument("-C", "--Acov", type=int, default=5, 
                    help='minimum A coverage for m6Am sites (default: 5)')
parser.add_argument("-s", "--Signal_Ratio", type=float, default=0.8, 
                    help="minimum ratio of signal reads (eg. reads with unconverted As less than 3) (default: 0.8)")
parser.add_argument("-r", "--methyl_Ratio", type = float, default=0.1, 
                    help="minimum m6A level (default: 0.1)")
parser.add_argument("-adp", "--FDR", type=float, default=0.05, 
                    help="FDR cutoff (default: 0.05)")
parser.add_argument("--persample", action="store_true", help="write out passed sites for each sample")
parser.add_argument("--skipStep1", action="store_true", 
                help="skip Step1, the merged sites list \{output\}.used must exists!")


args = parser.parse_args()
if args.prx is None:
    prx=[re.match("(.*/)?([^/]+).totalm6A.FDR.csv(.gz)?$", s).group(2) for s in args.inputs]
else:
    prx=args.prx
outdir = str(Path(args.output).parent)
outprx = str(Path(args.output).name)


# I. read file and filter sites with enough AG coverage, get a merged sites list --------
## 1. sites in the list with enough AG coverage in any sample
## 2. column "Passed" indicates a "ture" m6A site in any sample
## 3. "{output}.stat" file includes global m6A ratio for samples

if not args.skipStep1:
    print("[%s] Step1: Get reference sites list ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    for i in range(len(args.inputs)):
        print("\t[%s] Read sample %s ========" % 
            (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), prx[i]), 
            flush=True)

        # read file and filter sites
        df = pl.read_csv(args.inputs[i], separator="\t")
        if "Signal_Ratio" not in df.columns:
            df = df.with_columns(pl.lit(1).alias("Signal_Ratio"))
        df = df.filter(
            pl.col("AGcov") >= args.AGcov  # only sites with enough AG coverage used
        ).with_columns(
            ((pl.col("Acov") >= args.Acov) & (pl.col("Signal_Ratio") >= args.Signal_Ratio) & \
            (pl.col("Ratio") >= args.methyl_Ratio) & (pl.col("P_adjust") < args.FDR)).alias("Passed")
        )

        # write out passed sites for each sample
        if args.persample:
            df.filter(pl.col("Passed")).write_csv(outdir + "/" + prx[i] + "_" + outprx + ".m6A.passed", separator='\t')

        # global m6A ratio using A sites in each sample, respectively
        df = df.select(
            pl.col("Chr", "Sites", "Strand", "Gene", "AGcov", "Acov", "Passed")
        )
        df_stat = df.select(
            pl.col("AGcov").sum(),
            pl.col("Acov").filter(pl.col("Passed")).sum()
        ).with_columns(
            (pl.col("Acov")/pl.col("AGcov")*100).alias("Ratio"),
            pl.lit(prx[i]).alias("Sample")
        )

        # outer join samples
        if i == 0:
            df_merged = df.select(pl.col("Chr", "Sites", "Strand", "Gene", "Passed"))
            df_merged_stat = df_stat
        else:
            df_merged = df_merged.join(
                df.select(pl.col("Chr", "Sites", "Strand", "Gene", "Passed")), 
                on=["Chr", "Sites", "Strand", "Gene"], how="full", coalesce=True
            ).fill_null(False).with_columns(
                (pl.col("Passed") | pl.col("Passed_right")).alias("Passed")
            ).select(pl.col("Chr", "Sites", "Strand", "Gene", "Passed"))
            df_merged_stat = pl.concat([df_merged_stat, df_stat])

    # modify column name
    df_merged = df_merged.rename({"Sites": "Pos"})

    # output merged list
    df_merged.write_csv(args.output + ".m6A.used", separator='\t', null_value=".")
    df_merged_stat.write_csv(args.output + ".stat.sep", separator='\t')
else:
    df_merged = pl.read_csv(args.output + ".m6A.used", separator='\t')

# II. get passed sites with coverage info
print("\n[%s] Step2: Get merged passed sites info ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
df_passed = df_merged.filter(pl.col("Passed"))

for i in range(len(args.inputs)):
    print("\t[%s] Read sample %s ========" % 
          (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), prx[i]), 
          flush=True)

    df = pl.read_csv(args.inputs[i], separator="\t")
    if "Signal_Ratio" not in df.columns:
        df = df.with_columns(pl.lit(1).alias("Signal_Ratio"))
    
    # get used sites
    df = df_merged.join(
        df.rename({"Sites": "Pos"}), 
        on=["Chr", "Pos", "Strand", "Gene"], how="left", coalesce=True
    ).select(
        pl.col("Chr", "Pos", "Strand", "Gene", "Passed"),
        pl.col("AGcov").alias("AGcov_"+prx[i]),
        pl.col("Acov").alias("Acov_"+prx[i]),
        ((pl.col("AGcov") >= args.AGcov) & (pl.col("Acov") >= args.Acov) & (pl.col("Signal_Ratio") >= args.Signal_Ratio) & \
        (pl.col("Ratio") >= args.methyl_Ratio) & (pl.col("P_adjust") < args.FDR)).alias("Passed_"+prx[i])
    ).fill_null(False)

    # get global m6A ratio using common sites list in all samples
    df_stat2 = df.select(
        pl.col("AGcov_"+prx[i]).sum().alias("AGcov"),
        pl.col("Acov_"+prx[i]).filter(pl.col("Passed")).sum().alias("Acov")
    ).with_columns(
        (pl.col("Acov")/pl.col("AGcov")*100).alias("Ratio"),
        pl.lit(prx[i]).alias("Sample")
    )

    # filter common passed sites
    if i == 0:
        df_merged_stat2 = df_stat2
    else:
        df_merged_stat2 = pl.concat([df_merged_stat2, df_stat2])

    df_passed = df_passed.join(
        df,
        on=["Chr", "Pos", "Strand", "Gene", "Passed"], how="left", coalesce=True
    )


df_passed.write_csv(args.output + ".m6A.passed", separator="\t", null_value=".")
df_merged_stat2.write_csv(args.output + ".stat.common", separator='\t')

print("\n[%s] Done! ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)