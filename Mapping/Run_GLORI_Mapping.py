"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to run GLORI-tools"""
"""Input: [.fastq]"""

# version information
__version__ = "0.1.2"

"""
History:
    mod_v1.2: 
        1) add parameters to control pipeline
    mod_v1.1: 240326
        1) use awk to do A-G change
        2) make '-g' parameter usable for m6A_caller.py
    mod_v1.0: 
        1) Add comments; 
        2) not delete Important tmp files!!
        3) When coverge is < 15 for specifi segment (eg. chrM), pandas throw an error when output *.totalCR.txt, add an if-else to advoid.
        4) reformat log for timming
"""

import pandas as pd
import os, sys
import argparse
import time
import subprocess
import os
from time import strftime


def run_command(file,combine,untreated,rvs_fac,Threads):
    file3 = outputprefix + "_tf_rs.bam"  # mapped to transcriptome, G reversed to A, sorted
    file3_2 = outputprefix + "_rvs_rs.bam"  # mapped to reverse genome, G reversed to A, sorted
    file4 = outputprefix + ".trans2Genome.bam"  # transcriptome remapped to genome
    file5 = outputprefix + "_rs.bam"  # mapped to genome, G reversed to A, sorted
    file6 = outputprefix + ".trans2Genome.sorted.bam"  # transcriptome remapped to genome, sorted
    file6_2 = outputprefix + "_merged_t.bam"  # tmp merged
    file7_1 = outputprefix + "_merged.bam"  # merged
    finalbam = outputprefix + "_merged.sorted.bam"  # merged, sorted, final bam

    mapping_1 = "python "+DUOdir+"mapping_reads.py -i " + DUOdir + " -q "+ file +" -p "+ Threads + " -f "+ genome+ ' --FilterN '+FilterN
    mapping_2 = " -mulMax " + mulMax + " -t " + tool + " -m " + mismatch +" -pre "+ prx+ " -o " + outputdir
    if untreated:
        file3 = outputprefix + "_un_s.bam"
        file5 = outputprefix + "_s.bam"
        if combine:
            mapping_command = mapping_1 + " -Tf " + transgenome + mapping_2 + " --untreated --combine "
        else:
            mapping_command = mapping_1 + mapping_2 + " --untreated "
        if combine:
            print("\n**************combine,untreated", flush=True)
            print(mapping_command, flush=True)
            subprocess.call(mapping_command, shell=True)

            print("\n---- [%s] Transcript sites to genome locus " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            print("python " + DUOdir + "trans2genome.py --input " + file3 + " --output " + file4 + " --anno " + anno + " --fasta " + genome + " --sort --index", flush=True)
            subprocess.call("python " + DUOdir + "trans2genome.py --input " + file3 + " --output " + file4 +\
                            " --anno " + anno + " --fasta " + genome + " --untreated --sort --index",
                shell=True)
            print("\n---- [%s] Concat bam " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            subprocess.call("python " + DUOdir + "concat_bam.py -i " + file5 + " " + file6 + " -o " + file7_1 + " -t " + Threads + " --sort --index ",
                shell=True)         
            # subprocess.call("rm -f " + outputprefix +"_un*", shell=True)
            subprocess.call("rm -f " + outputprefix +".SJ.out.tab", shell=True)
            subprocess.call("rm -f " + outputprefix +".trans2Genome*", shell=True)
        else:
            print("\n**************uncombine,untreated", flush=True)
            print(mapping_command, flush=True)
            subprocess.call(mapping_command, shell=True)
            print("samtools view -F 4 -@ " + Threads + " -bS -h " + file5 + " | samtools sort -@ " + Threads + " > " + finalbam, flush=True)
            subprocess.call("samtools view -F 4 -@ " + Threads + " -bS -h " + file5 + " | samtools sort -@ " + Threads + " > " + finalbam,
                shell=True)
            subprocess.call("samtools index " + finalbam, shell=True)

        subprocess.call("rm -f " + outputprefix +"_s.bam*", shell=True)
        subprocess.call("mkdir -p "+outputdir+"/mapping-info", shell=True)
        subprocess.call("mv "+outputprefix+"*out "+outputdir+"/mapping-info", shell=True)
        subprocess.call("mv "+outputprefix+"*put "+outputdir+"/mapping-info", shell=True)
        sys.exit(0)
    else:
        if combine and rvs_fac:
            print("\n**************combine,treated", flush=True)
            mapping_command = mapping_1 + " -rvs " + rvsref +" -Tf "+ transgenome + mapping_2 + " --combine "+ " --rvs_fac"
            subprocess.call(mapping_command, shell=True)
            
            print("\n---- [%s] Transcript sites to genome locus " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            print("python "+DUOdir+"trans2genome.py --input " + file3 + " --output "+file4 + " --anno "+anno+" --fasta "+genome+" --sort --index", flush=True)
            subprocess.call("python "+DUOdir+"trans2genome.py --input " + file3 + " --output "+file4 + " --anno "+anno+" --fasta "+genome+" --sort --index",shell=True)
            print("\n---- [%s] Concat bam " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            subprocess.call("python "+DUOdir+"concat_bam.py -i " + file3_2 + " " + file5 + " -o " + file6_2 + " -t " + Threads + " --sort --index ",shell=True)
            filexx = file6_2[:-4]+'.sorted.bam'
            subprocess.call("python "+DUOdir+"concat_bam.py -i " + filexx + " " + file6 + " -o " + file7_1 + " -t " + Threads + " --sort --index ",shell=True)
            
            subprocess.call("rm -f " + filexx+"*", shell=True)
            subprocess.call("rm -f " + outputprefix +"*_un_2.fq", shell=True)
            subprocess.call("rm -f " + outputprefix +"*.SJ.out.tab", shell=True)
            subprocess.call("rm -f " + file3_2+"*", shell=True)
            subprocess.call("rm -f " + file3+"*", shell=True)
            subprocess.call("rm -f " + file5+"*", shell=True)
            subprocess.call("rm -f " + file6+"*", shell=True)
            subprocess.call("rm -f " + file7_1+"*", shell=True)
            subprocess.call("rm -f " + file6_2+"*", shell=True)
            subprocess.call("rm -f " + outputprefix + "_tf_rs.unlift.bam"+"*", shell=True)
        elif not combine and not rvs_fac:
            print("\n**************uncombine,treated,no rvs_fac", flush=True)
            mapping_command = mapping_1 + mapping_2
            print(mapping_command, flush=True)
            subprocess.call(mapping_command, shell=True)
            print("samtools view -F 20 -@ " + Threads + " -bS -h " + file5 + " | samtools sort -@ " + Threads + " > " + finalbam, flush=True)
            subprocess.call("samtools view -F 20 -@ " + Threads + " -bS -h " + file5 + " | samtools sort -@ " + Threads + " > " + finalbam,shell=True)
            subprocess.call("samtools index " + finalbam, shell=True)

        subprocess.call("mkdir -p "+outputdir+"/mapping-info", shell=True)
        subprocess.call("mv "+outputprefix+"*out "+outputdir+"/mapping-info", shell=True)
        subprocess.call("mv "+outputprefix+"*put "+outputdir+"/mapping-info", shell=True)
    print("\n---- [%s] merged bam finished " % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GLORI-tools for detecting m6A sites from high-throughput sequencing data")

    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-q", "--fastq", nargs="?", type=str, default=sys.stdin,
                        help="fastq files with surfix as _1.fq;_1.fastq;_2.fq;_2.fastq")
    group_required.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, help="Index file for the plus strand of the genome")
    group_required.add_argument("-f2", "--reference2", nargs="?", type=str, default=sys.stdin, help="Index file for the unchanged genome")
    group_required.add_argument("-rvs", "--rvsref", nargs="?", type=str, default=sys.stdin, help = "Index file for the minus strand of the genome")
    group_required.add_argument("-Tf", "--transref", nargs="?", type=str, default='None',help="Index file for the minus strand of the transcriptome")
    group_required.add_argument("-a", "--anno", nargs="?", type=str, default='None', help="Annotation file within exons")
    group_required.add_argument("--combine", "--combine", help="Whether mapping to transcriptome",action="store_true")
    group_required.add_argument("--untreated", "--untreated", help="If the input is untreated",action="store_true")
    group_required.add_argument("--rvs_fac", "--rvs_fac",help="Whether to map to the reverse strand of the transcriptome", action="store_true")

    group_output = parser.add_argument_group("Output (optional)")
    group_output.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default', help="Output file prefix")
    group_output.add_argument("-o", "--outputdir", nargs="?", type=str, default='./', help="Output directory")

    group_mappingfilter = parser.add_argument_group("Mapping conditions (optional)")
    group_mappingfilter.add_argument("-F", "--FilterN", nargs="?", type=str, default=0.5, help="The setting for the STAR parameter --outFilterScoreMinOverLread")
    group_mappingfilter.add_argument("-t", "--tools", nargs="?", type=str, default='STAR',  choices=['STAR', 'bowtie', 'bowtie2'],
                                     help="We recommend using STAR for genome alignment and Bowtie for transcriptome alignment.")
    group_mappingfilter.add_argument("-T", "--Threads", nargs="?", type=str, default='1',help="Used threads")
    group_mappingfilter.add_argument("-mulMax", "--mulMax", nargs="?", type=int, default=1, 
                                help="Suppress all alignments if > <int> exist (for bowtie and STAR); " + 
                                "set to -1 to allow unlimited multimapping (eg., for spike-in references with similar sequences, only for bowtie)")
    group_mappingfilter.add_argument("-m", "--mismatch", nargs="?", type=str, default='2', help="Permitted mapping mismatches")


    options = parser.parse_args()
    global genome,genome2,transgenome, outputdir,tool,Threads,mulMax,mismatch,prx,rvsref,outputprefix
    DUOdir = os.path.dirname(__file__)+"/"
    genome = options.reference
    genome2 = options.reference2
    transgenome = options.transref
    rvsref = options.rvsref
    outputdir = options.outputdir
    tool = options.tools
    Threads = options.Threads
    mulMax = options.mulMax
    mismatch = options.mismatch
    anno = options.anno
    prx = options.outname_prefix
    FilterN = str(options.FilterN)

    outputprefix = outputdir + "/" + prx
    if options.untreated:
        genome = genome.split(".AG_conversion.fa")[0]
        transgenome = transgenome.split(".AG_conversion.fa")[0]

    run_command(options.fastq,options.combine,options.untreated,options.rvs_fac,Threads)


