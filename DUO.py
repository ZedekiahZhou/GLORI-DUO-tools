#!/usr/bin/env python
"""
Author: Zhe Zhou, Peking University, Yi lab
Date: May 24, 2024
Email: zzhou24@pku.edu.cn
Program: This program is used for DUO-seq analysis
Version: 1.0.1
ToDo: 
    1. 最后一步call m6A可以考虑从拆分formatted_mpi开始， 而非referbase.mpi开始
    2. m6A样本的所有点信息？怎么和merge结合
"""

import argparse, subprocess, re, os, time, sys

def parse_unit_value(value_str):
    unit_map = {
        'K': 1000,
        'M': 1000000,
        'G': 1000000000,
    }

    if value_str[-1].upper() in unit_map:
        unit = value_str[-1].upper()
        number = float(value_str[:-1])
        return int(number * unit_map[unit])
    else:
        try:
            return int(value_str)
        except ValueError:
            raise argparse.ArgumentTypeError(f"Invalid argument '{value_str}': must be a numeric value with optional unit (K, M, G).")

parser = argparse.ArgumentParser(description="DUO-tools for detecting m6Am sites from DUO-seq data")

group_global = parser.add_argument_group("Global")
group_global.add_argument("--mode", type=str, default="m6Am", choices=["m6Am", "m6A"], 
                          help="run in m6Am or m6A mode, default is m6Am")
group_global.add_argument("--module", type=str, default=None, 
                          help="the module(s) to run: [preprocessing, mapping, call_m6Am, call_m6A or QC]; \
                            separated by comma [no space] if run multiple modules, default is running the whole pipeline")
group_global.add_argument("--untreated", action="store_true", 
                          help="untreated samaples (only preprocessing and mapping)")
group_global.add_argument("--raw_fq", type=str, help="raw fastq(.gz) file")
group_global.add_argument("--clean_fq", type=str, help="clean fastq file for mapping")
group_global.add_argument("--bam", type=str, help="merged sorted bam files")
group_global.add_argument("--prx", type=str, help="output file prefix")
group_global.add_argument("--DUOdir", type=str, default=os.path.dirname(__file__)+"/", help="directory of DUO-tools")
group_global.add_argument("-o", "--outdir", type=str, default="./",
                          help="Default is current directory. \
                            Outputs of preprocessing will be saved in {outdir}/02_Clean/; \
                            outputs of fastqc will be saved in {outdir}/fastq/; \
                            outputs of mapping and calling will be saved in {outdir}/03_Sites/; \
                            and outputs of QC will be saved in {outdir}/03_Sites/QC/.")
group_global.add_argument("-p", "--threads", type=int, default=20, help="threads used, default is 20")
group_global.add_argument("--test", action="store_true", help="only print commands to run")

group_preprocessing = parser.add_argument_group("Preprocessing")
group_preprocessing.add_argument("--umi5", dest="umi5", type=int, default=0, 
                                 help="length of bases to removed from 5 prime of the reads, default is 0")
group_preprocessing.add_argument("--umi3", dest="umi3", type=int, default=0, 
                                 help="length of bases to removed from 3 prime of the reads, default is 0")
group_preprocessing.add_argument("--nextseq", action="store_true", 
                                 help="NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.")
group_preprocessing.add_argument("-m", "--min_len", type=int, default=25, 
                                 help="discard reads that bacame shorter than min_len before mapping")
group_preprocessing.add_argument("--fastqc", dest="fastqc", type=bool, default=True, 
                                 help="whether to run fastqc, default is True")
group_preprocessing.add_argument("--tag_seq", dest="tag_seq", type=str, default="TGACGCTGCCGACGATC", 
                                 help="tag sequence ligated to 5' ends of TSS, default is TGACGCTGCCGACGATC. \
                                    Skip tag removal if set to 'none'.")

group_mapping = parser.add_argument_group("Mapping")
group_mapping.add_argument("-f", "--reference", nargs="?", help="Index file for the plus strand of the genome")
group_mapping.add_argument("-f2", "--reference2", nargs="*", help="Index file for the unchanged genome")
group_mapping.add_argument("-rvs", "--rvsref", nargs="?", help = "Index file for the minus strand of the genome")
group_mapping.add_argument("-Tf", "--transref", nargs="?", help="Index file for the minus strand of the transcriptome")
group_mapping.add_argument("-a", "--anno", nargs="?", help="Annotation file within exons")
    
group_m6Am = parser.add_argument_group("Call m6Am or m6A")
group_m6Am.add_argument("-ds", "--ds2N", type=parse_unit_value, default=None, 
                        help="[Both] Downsample the merged sorted bam files to N reads, default is no downsample.")

group_m6Am.add_argument("-c", "--cov", type=int, default=15, 
                        help='[Both] minimum A+G coverage for TSS or m6A sites, default is 15')
group_m6Am.add_argument("-C", "--Acov", type=int, default=5, 
                        help='[Both] minimum A coverage for m6A(m) sites, default is 5')
group_m6Am.add_argument("-s", "--Signal_Ratio", type=float, default=0.8, 
                        help="[Both] minimum ratio of signal reads (eg. reads with unconverted As less than 3), \
                            default is 0.8")
group_m6Am.add_argument("-R", "--AG_Ratio", type=float, default=0.8, 
                    help="[Both] minimum ratio of (A+G reads)/total in this sites, default is 0.8")
group_m6Am.add_argument("-adp", "--FDR", type=float, default=0.001, help="[Both] FDR cutoff, default is 0.001")
group_m6Am.add_argument("-ta", "--tssanno", nargs="*", type=str, help="[m6Am] Annotation of TSS range")
group_m6Am.add_argument("--tpm", type=float, default=1.0, help="[m6Am] minimum TPM value for TSS, default is 1.0")
group_m6Am.add_argument("--absDist", type=int, default=1000, 
                        help="[m6Am] maximum absolute distance to any annotated TSS from GTF file, default is 1000")
group_m6Am.add_argument("--prop", type=float, default=0.05,
                        help="[m6Am] minimum proportion relative to the total TPM of a gene")
group_m6Am.add_argument("--zscore", type=float, default=1.0,
                        help="[m6Am] minimum Z-score (calculated within a gene) for TSS, default is 1.0")
group_m6Am.add_argument("-ba", "--baseanno", type=str, default='None',help="[m6A] Annotations at single-base resolution")
group_m6Am.add_argument("-r", "--methyl_Ratio", type = float, default=0.1, help="[m6A] minimum m6A level")

group_QC = parser.add_argument_group("QC")
group_QC.add_argument("--gtf", type=str, help="GTF file, chromosome name with suffix '_AG_converted'")
group_QC.add_argument("--gtf2", type=str, help="GTF file")
group_QC.add_argument("--mPRdir", type=str, help = "Directory of metaPlotR")
group_QC.add_argument("--mPRanno", type=str, help="Directory and prefix of metaPlotR annotations")

args = parser.parse_args()

def run_cmd(cmd):
    global args
    print("--- [%s] " % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + cmd, flush=True)
    if not args.test: subprocess.call(cmd, shell=True)


def fun_pre(raw_fq, prx, args):
    """
    Preprocessing raw fastq file (remove adapter, duplicates and 5' tag)
    """
    print(prx, flush=True)
    print("\n[%s] Preprocessing ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)

    # fastqc on raw fastq
    if args.fastqc:
        if not args.test: 
            subprocess.call("mkdir -p " + args.outdir +"/fastqc/raw " + args.outdir + "/fastqc/clean", shell=True)
        run_cmd("fastqc "+ raw_fq + " --thread " + str(args.threads) + " -q -o " + args.outdir +"/fastqc/raw")

    clean_dir=args.outdir + "/02_Clean/"

    # trim adapter
    cutoff_len1=args.umi5+args.umi3
    if args.nextseq:
        run_cmd("trim_galore --nextseq 20 -j 7 --stringency 1 -e 0.3 --length " + str(cutoff_len1) + 
                " -o " + clean_dir + " " + raw_fq)
    else:
        run_cmd("trim_galore -q 20 -j 7 --stringency 1 -e 0.3 --length " + str(cutoff_len1) + 
                " -o " + clean_dir + " " + raw_fq)

    tmpfile=clean_dir + re.match(".+/([^/]+).f(ast)?q(.gz)?$", raw_fq).group(1) + "_trimmed.fq.gz"
    trimmed_fq=clean_dir + prx + "_trimmed.fq.gz"
    run_cmd("mv " + tmpfile + " " + trimmed_fq)

    # remove duplicates and umi
    rmdup_fq=clean_dir + prx + "_rmdup.fq.gz"
    run_cmd("seqkit rmdup -j 10 -s " + trimmed_fq + " | fastx_trimmer -Q 33 -f " + str(args.umi5+1) + \
            " -z -o " + rmdup_fq)

    # remove tag
    clean_fq=clean_dir + prx + "_clean.fq"
    tag_info=clean_dir + prx + "_rmtag.info"
    if args.tag_seq.lower() == "none":
        # only remove reads shorter than min_len
        cmd='cutadapt -j 0 -m ' + str(args.min_len) + ' -o ' + clean_fq + ' ' + rmdup_fq
        run_cmd(cmd)
    else:
        cmd='cutadapt -j 0 -g "' + args.tag_seq + ';rightmost" -m ' + str(args.min_len) + \
            ' -O ' + str(len(args.tag_seq)) + ' -e 0.2 -o ' + clean_fq + ' --info-file ' + tag_info + ' ' + rmdup_fq
        if args.mode == "m6Am":
            cmd = cmd + ' --discard-untrimmed'
        run_cmd(cmd)
    
    # fastqc on clean fastq
    if args.fastqc:
        run_cmd("fastqc "+ clean_fq + " --thread " + str(args.threads) + " -q -o " + args.outdir +"/fastqc/clean")


def fun_mapping(clean_fq, prx, args):
    """
    Ternary mapping using STAR
    """
    print(prx, flush=True)
    print("\n[%s] Mapping ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)

    site_dir=args.outdir + "/03_Sites/"
    if args.reference is None:
        raise ValueError("Please provide the plus strand reference (-f) file for mapping.")
    if args.transref is None:
        raise ValueError("Please provide the transcriptome reference (-Tf) file for mapping.")
    cmd_prx="python " + args.DUOdir + "/Mapping/Run_GLORI_Mapping.py -q " + clean_fq + " -T " + str(args.threads) + \
        " -f " + args.reference + " -pre " + prx + " -o " + site_dir
    if args.untreated:
        cmd=cmd_prx + " -Tf " + args.transref + " -a " + args.anno + " --combine --untreated"
        run_cmd(cmd)
    else:
        if args.reference2 is None:
            raise ValueError("Please provide the unchanged reference (-f2) file for mapping GLORI treated samples.")
        if args.rvsref is None:
            raise ValueError("Please provide the minus strand reference (-rvs) file for mapping GLORI treated samples.")
        cmd=cmd_prx + " -f2 " + " ".join(args.reference2) + " -rvs " + args.rvsref+ " -Tf " + args.transref + \
            " -a " + args.anno + " --combine --rvs_fac"
        run_cmd(cmd)
        run_cmd("mv " + site_dir + prx + "_A.bed_sorted " + site_dir + prx + "_AGchanged_2.fq " + \
                site_dir + "/intermediate/")

    # fastqc on mapped bam files
    if args.fastqc:
        if not args.test: 
            subprocess.call("mkdir -p " + args.outdir +"/fastqc/mapped ", shell=True)
        mapped_bam=args.outdir + "/03_Sites/" + prx + "_merged.sorted.bam"
        run_cmd("fastqc " + mapped_bam + " --thread " + str(args.threads) + " -q -o " + args.outdir +"/fastqc/mapped")
        


def fun_m6Am(bam, prx, args):
    """
    Call m6Am
    """
    print(prx, flush=True)
    print("\n[%s] m6am calling ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    site_dir=args.outdir + "/03_Sites/"

    # downsample
    if args.ds2N is not None:
        n_reads=subprocess.run("samtools view -c " + bam, shell=True, capture_output=True, text=True)
        dsr=int(args.ds2N)/int(n_reads.stdout)
        if dsr>=1: dsr=0.99999

        prx=prx+"_ds"
        ds_bam = site_dir + prx + "_merged.sorted.bam"
        run_cmd("samtools view -@ 4 -s " + str(dsr) + " " + bam + " -b > " + ds_bam + " && samtools index " + ds_bam)

        bam=ds_bam  # use ds_bam for downstream analysis

    # call TSS
    ftss=site_dir + prx + "_TSS_raw.bed"
    run_cmd("python " + args.DUOdir + "/Call_m6Am/get_TSS.py -r " + " ".join(args.reference2) + " -b " + bam + " -o " + ftss)

    # annotate TSS
    ftss_anno=site_dir + prx + "_TSS_raw.bed.annotated"
    run_cmd("bedtools intersect -a " + ftss + " -b " + " ".join(args.tssanno) + " -s -wa -wb -loj > " + ftss_anno)
    run_cmd("python " + args.DUOdir + "/Call_m6Am/anno_TSS.py -i " + ftss_anno) # output to {ftss_anno}.rmdup

    if not args.untreated:
        # pile ATCG
        fAGcount=site_dir + prx + "_AGcount.tsv"
        run_cmd("python " + args.DUOdir + "/Call_m6Am/pileup_reads5p_multiprocessing.py -p " + str(args.threads) + \
                " -r " + " ".join(args.reference2) + \
                " -l " + ftss_anno + ".rmdup -o " + fAGcount + " -b " + bam)

        # call m6Am
        run_cmd("python " + args.DUOdir + "/Call_m6Am/m6Am_caller.py -i " + fAGcount)
        
        #run_cmd("mv " + fAGcount + " " + site_dir + prx + "_m6Am_sites_raw.tsv " + site_dir + "/intermediate/")
    
    run_cmd("rm -rf " + ftss_anno)
    run_cmd("mv " + ftss + " " + site_dir + "/intermediate/")
    #run_cmd("mv " + ftss_anno + ".rmdup " + site_dir + "/intermediate/")
    

def fun_m6A(bam, prx, args):
    """
    Call m6A
    """
    print(prx, flush=True)
    print("\n[%s] m6a calling ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    site_dir=args.outdir + "/03_Sites/"

    # downsample
    if args.ds2N is not None:
        n_reads=subprocess.run("samtools view -c " + bam, shell=True, capture_output=True, text=True)
        dsr=int(args.ds2N)/int(n_reads.stdout)
        if dsr>=1: dsr=0.99999

        prx=prx+"_ds"
        ds_bam = site_dir + prx + "_merged.sorted.bam"
        run_cmd("samtools view -@ 4 -s " + str(dsr) + " " + bam + " -b > " + ds_bam + " && samtools index " + ds_bam)

        bam=ds_bam  # use ds_bam for downstream analysis
    
    # pileup
    run_cmd("python " + args.DUOdir + "/Call_m6A/Run_GLORI_pileup.py --bam " + bam + " -T " + str(args.threads) + \
        " -f " + args.reference + " -f2 " + " ".join(args.reference2) + " -pre " + prx + " -o " + site_dir)

    # call m6A
    fmpi=site_dir + prx + ".referbase.mpi"
    run_cmd("python " + args.DUOdir + "/Call_m6A/Run_GLORI_m6A.py --mpi " + fmpi + " -T " + str(args.threads) + \
        " -b " + args.baseanno + " -c " + str(args.cov) + " -C " + str(args.Acov) + " -adp " + str(args.FDR) + \
        " -s " + str(args.Signal_Ratio) + " -r " + str(args.methyl_Ratio) + " -R " + str(args.AG_Ratio) + \
        " -pre " + prx + " -o " + site_dir)


def fun_QC(bam, prx, args):
    """
    Quality control
    """
    print(prx, flush=True)
    print("\n[%s] Quality control ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
    qc_dir=args.outdir + "/03_Sites/QC/"

    # 1. get alignments length
    print("--- [%s] " % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + "Get alignments length:", flush=True)
    mini_bam = qc_dir + prx + "_subsample_0.05.bam"
    run_cmd("samtools view -@ 4 -s 0.05 " + bam + " -b > " + mini_bam + " && samtools index " + mini_bam)
    run_cmd("python " + args.DUOdir + "/QC/get_align_len.py " + mini_bam + " " + qc_dir + prx + "_align_len.csv") 

    # 2. get read distribution (qualimap)
    print("--- [%s] " % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + "Get read distribution (qualimap):", flush=True)
    qualimap_dir = qc_dir + "qualimap/" + prx + "/"
    gtf = args.gtf if not args.untreated else args.gtf2
    run_cmd("mkdir -p " + qualimap_dir + " && qualimap rnaseq -outdir " + qualimap_dir + " -bam " + mini_bam + " -gtf " + gtf + " --java-mem-size=32G")


def main(args):
    os.makedirs(args.outdir+"/02_Clean/", exist_ok=True)
    os.makedirs(args.outdir+"/03_Sites/QC/", exist_ok=True)
    os.makedirs(args.outdir+"/03_Sites/intermediate/", exist_ok=True)

    # parse modules to run with
    if args.module is None:
        module = ["preprocessing", "mapping", "QC"]
        if args.mode == "m6Am":
            module += ["call_m6Am"]
        elif not args.untreated:
            module += ["call_m6A"]
    else:
        module = args.module.split(",")
        if args.mode == "m6Am":
            if "call_m6A" in module:
                raise ValueError("No call_m6A module in m6Am mode!")
        else:
            if "call_m6Am" in module:
                raise ValueError("No call_m6Am module in m6A mode!")
    
    # call each module
    prx=None  # Initialize prx to detect the beginning module

    # Check:
    #   1. if the pipelines are begin with spcific steps
    #   2. wheather args.prx were specified by user

    if "preprocessing" in module:
        if args.raw_fq is None:
            raise ValueError("Raw fastq files must be provided if run preprocessing!")
        if args.prx is None:
            prx=re.match("(.*/)?([^/]+)_[LS][0-9]+_(R)?[12].f(ast)?q(.gz)?$", args.raw_fq).group(2)
        else:
            prx=args.prx
        fun_pre(args.raw_fq, prx, args)

    if "mapping" in module:
        if prx is None:  # begin with mapping step
            if args.clean_fq is None:
                raise ValueError("Clean fastq files must be provided if beginning with the mapping step!")
            else:
                if args.prx is None:
                    prx=re.match("(.*/)?([^/]+)_clean.fq$", args.clean_fq).group(2)
                else:
                    prx=args.prx
            clean_fq=args.clean_fq
        else:
            clean_fq=args.outdir + "/02_Clean/" + prx + "_clean.fq"
        fun_mapping(clean_fq, prx, args)

    if "call_m6Am" in module:
        if prx is None:  # begin with m6Am calling step
            if args.bam is None:
                raise ValueError("Bam files must be provided if beginning with the m6Am calling step!")
            else:
                if args.prx is None:
                    prx=re.match("(.*/)?([^/]+)_merged.sorted.bam$", args.bam).group(2)
                else:
                    prx=args.prx
            bam=args.bam
        else:
            bam=args.outdir + "/03_Sites/" + prx + "_merged.sorted.bam"
        fun_m6Am(bam, prx, args)

    if "call_m6A" in module:
        if prx is None:  # begin with m6A calling step
            if args.bam is None:
                raise ValueError("Bam files must be provided if beginning with the m6A calling step!")
            else:
                if args.prx is None:
                    prx=re.match("(.*/)?([^/]+)_merged.sorted.bam$", args.bam).group(2)
                else:
                    prx=args.prx
            bam=args.bam
        else:
            bam=args.outdir + "/03_Sites/" + prx + "_merged.sorted.bam"
        fun_m6A(bam, prx, args)
    
    if "QC" in module:
        try:
            bam
        except NameError:
            bam=None
        if prx is None:  # begin with QC step
            if args.bam is None:
                raise ValueError("Bam files must be provided if beginning with the QC step!")
            else:
                if args.prx is None:
                    prx=re.match("(.*/)?([^/]+)_merged.sorted.bam$", args.bam).group(2)
                else:
                    prx=args.prx
            if bam is None:
                bam=args.bam
        else:
            if bam is None:
                bam=args.outdir + "/03_Sites/" + prx + "_merged.sorted.bam"
        fun_QC(bam, prx, args)

    print("\n[%s] Done! ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


if __name__ == '__main__':
    main(args)
