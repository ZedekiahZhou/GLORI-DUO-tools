import polars as pl
from scipy import stats
import argparse
from statsmodels.stats.multitest import multipletests

"""
History:
    1. remove "-p", do not filter p value, use all rows to calculate FDR in "m6A_caller_FDRfilter.py"
"""

# args
parser = argparse.ArgumentParser(description="call m6Am sites from AGcounts file")
parser.add_argument("-i", "--input", type=str, required=True, help="AGcounts file from pileup_reads5p.py")
parser.add_argument("-o", "--output", type=str, required=False, help="output file name")
args = parser.parse_args()

# input and output
df=pl.read_csv(args.input, separator="\t")

if args.output is None:
    args.output = args.input.replace("_AGcount.tsv", "") + "_m6Am_sites_raw.tsv"


# binomtest function
def test_sites(x):
    if x["AG_cov"] == 0:
        return 1
    else:
        return stats.binomtest(x["A"], n=x["AG_cov"], p=x["Ctrl_Ratio"], alternative="greater").pvalue


# perform test
df = df.with_columns(
    (pl.col("A") + pl.col("T") + pl.col("C") + pl.col("G")).alias("Signal_cov"),  # signal reads (with unconverted As <= 3?)
    (pl.col("A") + pl.col("G")).alias("AG_cov"),  # A + G coverage (signal)
    (pl.col("Next_pos_A") + pl.col("Next_pos_G")).alias("Next_pos_AG")  # Ctrl A + G coverage
).with_columns(
    (pl.col("Signal_cov")/pl.col("Counts")).round(3).alias("Signal_Ratio"),  # signal ratio
    (pl.col("AG_cov")/pl.col("Signal_cov")).round(3).alias("AG_Ratio"),  # AG ratio
    pl.when(pl.col("Next_pos_AG") == 0)
    .then(pl.lit(0))
    .otherwise(pl.col("Next_pos_A")/pl.col("Next_pos_AG"))
    .round(3).alias("Ctrl_Ratio"),  # Ctrol unconverted ratio
    pl.when(pl.col("AG_cov") == 0)
    .then(pl.lit(0))
    .otherwise(pl.col("A")/pl.col("AG_cov"))
    .round(3).alias("m6Am_Ratio")  # m6Am ratio
).with_columns(
    pl.struct(["A", "AG_cov", "Ctrl_Ratio"])
    .map_elements(test_sites, return_dtype=pl.Float64)
    .alias("Pvalue")  # binom_test p value
).with_columns(
    pl.col("Pvalue")
    .map_batches(lambda x: multipletests(x, method="fdr_bh")[1])
    .alias("FDR")  # BH adjusted
)

# columns order inherited from input, add columns Signal_cov, ...
df.write_csv(args.output, separator="\t")