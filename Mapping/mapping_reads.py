"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: Align reads to Genome and trancriptome. This script is modified from GLORI"""
"""Input: [fastq files]"""

# version information
__version__ = "0.1.2"

"""
History: Version evolved with main pipeline
    mod_v1.2:
        1) add bowtie2 for target-seq
    mod_v1.1:
        1) use awk to do A-G change
    mod_v1.0:
        1) reformat log for timming
"""

# import mxnet
import sys
import os
import re
import argparse
import subprocess
import itertools
import time
from heapq import merge
import glob
from time import strftime
from Bio.Seq import reverse_complement


parser = argparse.ArgumentParser(description = "reads alignment")
parser.add_argument("-i", "--DUOdir", nargs="?", type=str, default=sys.stdin,help = "Directory where DUO-tool is located")
parser.add_argument("-q", "--fastq", nargs="?", type=str, default=sys.stdin, help = "fastqfiles with surfix as _1.fq;_1.fastq;_2.fq;_2.fastq")
parser.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, help = "Index file for the genome")
parser.add_argument("-rvs", "--rvsref", nargs="?", type=str, default=sys.stdin, help = "Index file for the minus strand of the genome")
parser.add_argument("-Tf", "--transref", nargs="?", type=str, default=sys.stdin, help = "Index file for the minus strand of the transcriptome")
parser.add_argument("-t", "--tools", nargs="?", type=str, default=sys.stdin,
                    help="We recommend using STAR for genome alignment and Bowtie for transcriptome alignment")
parser.add_argument("-m", "--mismatch", nargs="?", type=int, default=2, help="Permitted mapping mismatches")
parser.add_argument("-F", "--FilterN", nargs="?", type=str, default=0.5, help="The setting for the STAR parameter --outFilterScoreMinOverLread")
parser.add_argument("-mulMax", "--mulMax", nargs="?", type=int, default=1, 
                    help="Suppress all alignments if > <int> exist (for bowtie and STAR); " + 
                    "set to -1 to allow unlimited multimapping (eg., for spike-in references with similar sequences, only for bowtie)")
parser.add_argument("--combine", "--combine", help="Whether mapping to transcriptome",action="store_true")
parser.add_argument("--untreated", "--untreated", help="If the input is untreated",action="store_true")
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',help = "--outname_prefix")
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default=sys.stdin, help="outputdir")
parser.add_argument("--rvs_fac", "--rvs_fac",help="Whether to map to the reverse strand of the transcriptome", action="store_true")
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "Used threads")
args = parser.parse_args()

DUOdir = args.DUOdir
fastq = args.fastq
Threads = args.Threads
reference = args.reference
rvsref = args.rvsref
transref = args.transref
tools = args.tools
global FilterN
FilterN = args.FilterN

mismatch = args.mismatch
mulMax = args.mulMax
outputdir = args.outputdir
outname_prx = args.outname_prefix


re_digits = re.compile(r'(\d+)')

def run_cmd(cmd):
    print("--- [%s] " % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + cmd, flush=True)
    subprocess.call(cmd, shell=True)

def embedded_numbers(s):
    s2=s.strip().split("\t")
    pieces = re_digits.split(s2[0])
    pieces[1::2] = map(int, pieces[1::2])
    return pieces


def sort_bedfiles(bedfiles,outputfiles):
    prx2 = bedfiles[:-4]
    path = prx2+"_chunk_*.bed"
    chunksize = 5000000
    fid = 1
    lines = []
    with open(bedfiles, 'r') as f_in:
        f_out = open(prx2+'_chunk_{}.bed'.format(fid), 'w')
        for line_num, line in enumerate(f_in, 1):
            lines.append(line)
            if not line_num % chunksize:
                lines = sorted(lines, key=embedded_numbers)
                f_out.writelines(lines)
                f_out.close()
                lines = []
                fid += 1
                f_out = open(prx2+'_chunk_{}.bed'.format(fid), 'w')
        # last chunk
        if lines:
            lines = sorted(lines, key=embedded_numbers)
            f_out.writelines(lines)
            f_out.close()
            lines = []

    chunks = []
    for filename in glob.glob(path):
        chunks += [open(filename, 'r')]

    with open(outputfiles, 'w') as f_out:
        print('merging bedfiles')
        f_out.writelines(merge(*chunks, key=embedded_numbers))
    subprocess.call("rm -f " + prx2 + "_chunk_*.bed", shell=True)


def mapping_files(tool,fastq,reference,Threads,muta_N,fqname,outputdir,mulMax,flag):
    outputfile = outputdir +fqname+".sam"
    unmapfastq = outputdir +fqname+"_un_2.fq"
    if tool == "bowtie2":
        command = "bowtie2 --local --no-unal -p " + Threads + " -x "+reference + " -U "+fastq + " -S "+outputfile + " --un " + unmapfastq + " 2>" + outputfile + '.output'
        run_cmd(command)
        run_cmd("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile + 
                " | samtools sort -n -O SAM > " + outputfile + "_sorted.sam")
        run_cmd("mv -f " + outputfile + "_sorted.sam " + outputfile)
    elif tool == "bowtie":
        if mulMax == -1:
            para_0 = 'bowtie -k 1'
        else:
            para_0 = 'bowtie -k 1 -m '+ str(mulMax)
        para_A = ' -v '+ str(muta_N)
        para_B = ' --best --strata -p ' + Threads
        para_C = ' -x '+ reference +" "+ fastq +' -S ' + outputfile
        para_unmap = ' --un ' + unmapfastq
        para_end = ' 2>' + outputfile +'.output'
        command = para_0+para_A+para_B+para_C+para_unmap+para_end
        run_cmd(command)
        run_cmd("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile +
                " | samtools sort -n -O SAM > " + outputfile + "_sorted.sam")
        run_cmd("mv -f " + outputfile + "_sorted.sam " + outputfile)
    elif tool == "STAR":
        para_0 = "STAR --runThreadN "+ Threads + " --readNameSeparator ' '"
        para_g = " --genomeDir "+ reference[:-3]
        para_A = " --limitOutSJcollapsed 5000000 "
        para_B = " --outFilterMismatchNmax " + str(muta_N)
        # para_B_2 = " --outFilterMismatchNoverLmax 0.3"
        # para_B_3 = " --outFilterMismatchNoverReadLmax 1"
        para_B_2=''
        para_B_3=' --outFilterScoreMinOverLread '+FilterN+' --outFilterMatchNminOverLread '+FilterN+' --seedSearchStartLmax 30 '# increase overall mapping sensitivity
        para_C = " --outSAMattributes All --outSAMprimaryFlag AllBestScore --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMtype BAM Unsorted"
        para_D = " --outFilterMultimapNmax " + str(mulMax)
        para_E = " --outFileNamePrefix " + outputfile[:-3] + " --readFilesIn " + fastq
        para_unmap = " --outSAMunmapped Within --outReadsUnmapped Fastx"
        line_command = para_0+para_g+para_A+para_B+para_B_2+para_B_3+para_C+para_D+para_E + para_unmap
        run_cmd(line_command)
        run_cmd("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile[:-3] + "Aligned.out.bam " + 
                "| samtools sort -n -O SAM > " + outputfile)
        run_cmd("mv " + outputfile[:-3] + 'Unmapped.out.mate1 ' + unmapfastq)
        run_cmd("rm -f " + outputfile[:-3] + 'Aligned.out.sam')
        run_cmd("rm -f " + outputfile[:-3] + 'Aligned.out.bam')
    return outputfile, unmapfastq


def getbamfiles(outputfile,fac,Threads,flag):
    output_bam = outputfile[:-4] + fac
    run_cmd("samtools view -F " + flag + " -bS -@ " + Threads+" -h "+outputfile + 
            " | samtools sort > " + output_bam)
    run_cmd("samtools index " + output_bam)
    run_cmd("rm -f " + outputfile)
    return output_bam


if __name__ == "__main__":
    global step
    step = 10000
    global change_fac,fqname2
    change_fac = 'AG'
    if outname_prx != 'default':
        fqname = outname_prx
    else:
        fqname = "_".join(os.path.basename(fastq).split(".")[:-1])
    outputdir2 = outputdir+"/"
    if os.path.exists(outputdir2):
        pass
    else:
        os.makedirs(outputdir2)

    fqname2= outname_prx

    if args.untreated:
        sys.stderr.write("[%s]untreated...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        print("\n---- [%s] Mapping to genome " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
        outputfile_untreated, unmapfastq = mapping_files(tools, fastq, reference, Threads, mismatch,
                                                        fqname2, outputdir2, mulMax,'4')
        untreated_bam = getbamfiles(outputfile_untreated,"_s.bam",Threads,'4')
        if args.combine:
            sys.stderr.write("[%s]untreated map to transcriptome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            print("\n---- [%s] Mapping to transcriptome " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            outputfile_untreated_unmap, _, = mapping_files('bowtie', unmapfastq,  transref, Threads,
                                                mismatch,fqname2 + "_un", outputdir2, mulMax,'4')
            bamAG_unmap = getbamfiles(outputfile_untreated_unmap,"_s.bam", Threads,'4')
    else:
        changefastq = outputdir2 + "/" + fqname2 + "_" + change_fac + "changed_2.fq"
        output_bed = outputdir2 + "/" + fqname2 + "_A.bed"
        if os.path.exists(changefastq):
            sys.stderr.write("[%s] Warning: the changed files already exists, please make sure the input file is correct\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            pass
        else:
            sys.stderr.write("[%s] change to A>G...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            print("\n---- [%s] Change A to G " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            subprocess.call("rm -f " + changefastq + " " + output_bed + " 2>/dev/null",shell=True)
            cmd = 'awk -f ' + DUOdir + 'AtoG.awk changefastq="' + changefastq + '" ' + fastq + ' > ' + output_bed
            run_cmd(cmd)

            print("---- [%s] Sort bed " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            sort_bedfiles(output_bed, output_bed + "_sorted")
            run_cmd("rm -f " + output_bed)

        sys.stderr.write("[%s] map to genome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        print("\n---- [%s] Mapping to genome " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
        outputfile_changeAG, unmapfastq = mapping_files(tools, changefastq, reference, Threads, mismatch,
                                                      fqname2, outputdir2,mulMax,'20')
        print("---- [%s] Reverse G to A " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
        reverse_samAG = outputfile_changeAG[:-4] + "_r.sam"
        cmd='awk -f '+ DUOdir + 'reverseGtoA.awk fbed="' + output_bed+'_sorted' + '" ' + outputfile_changeAG + ' > ' + reverse_samAG
        run_cmd(cmd)
        run_cmd("rm -f " + outputfile_changeAG)
        reversed_bamAG = getbamfiles(reverse_samAG,'s.bam', Threads,'20')

        if args.rvs_fac:
            sys.stderr.write("[%s]map to reversegenome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            print("\n---- [%s] Mapping to reversegenome " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            rvsmapfastq = outputdir2 + fqname2 + "_un_2.fq"
            outputfile_changeAG_rvsmap, _, = mapping_files(tools, rvsmapfastq, rvsref, Threads, mismatch,
                                                          fqname2 + "_rvs", outputdir2, mulMax, '4')
            print("---- [%s] Reverse G to A " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            reverse_samAG_rvsmap = outputfile_changeAG_rvsmap[:-4] + "_r.sam"
            cmd='awk -f '+ DUOdir + 'reverseGtoA_rvs.awk fbed="' + output_bed+'_sorted' + '" ' + outputfile_changeAG_rvsmap + ' > ' + reverse_samAG_rvsmap
            run_cmd(cmd)
            run_cmd("rm -f " + outputfile_changeAG_rvsmap)
            reversed_bamAG_rvsmap = getbamfiles(reverse_samAG_rvsmap, 's.bam', Threads, '4')
        else:
            pass

        if args.combine:
            unmapfastq = outputdir2 + fqname2 + "_rvs_un_2.fq"
            sys.stderr.write("[%s]map to transcriptome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            print("\n---- [%s] Mapping to transcriptome " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            outputfile_changeAG_unmap2, _, = mapping_files('bowtie', unmapfastq, transref, Threads, mismatch,
                                                          fqname2 + "_tf", outputdir2, mulMax,'20')
            print("\n---- [%s] Reverse G to A " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            reverse_samAG_unmap2 = outputfile_changeAG_unmap2[:-4] + "_r.sam"
            cmd='awk -f '+ DUOdir + 'reverseGtoA.awk fbed="' + output_bed+'_sorted' + '" ' + outputfile_changeAG_unmap2 + ' > ' + reverse_samAG_unmap2
            run_cmd(cmd)
            run_cmd("rm -f " + outputfile_changeAG_unmap2)
            reversed_bamAG_unmap2 = getbamfiles(reverse_samAG_unmap2, 's.bam', Threads,'20')