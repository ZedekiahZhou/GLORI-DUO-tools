# GLORI-DUO-tools

GLORI-DUO-tools is a bioinformatics toolkit designed for the analysis of GLORI-DUO sequencing data, which enables parallel quantification of m6A and m6Am modifications in the transcriptome.

The pipeline extends the original [GLORI-tools](https://github.com/liucongcas/GLORI-tools) framework and provides a one-command workflow for data preprocessing, read mapping, and modification site identification.

---

# Workflow

GLORI-DUO-tools includes two modes (`m6A` and `m6Am`) and four modules (`preprocessing`, `mapping`, `call_m6A[m]`, and `QC`) for each mode. These can be specified using the `--mode` and `--module` parameters.

## Important Notes

GLORI-DUO-tools currently supports single-end sequencing reads only.

Before running the pipeline, users should ensure that:

* reads are strand-specific and correctly oriented
* A-to-G conversion has been properly introduced

---

## System Requirements

To install GLORI-DUO-tools, simply download the source code and extract the files.

The following software packages are required:

| Software | Recommended Version |
| -------- | ------------------- |
| Python   | >= 3.8              |
| STAR     | >= 2.7              |
| bowtie   | >= 1.3              |
| samtools | >= 1.10             |
| bedtools | >= v2.31.1          |

The following Python packages are also required:

```text
polars
pysam
pandas
numpy
scipy
statsmodels
biopython
argparse
multiprocessing
collections
subprocess
glob
itertools
```

---

# Reference and Annotation Preparation

To prepare references and annotations for the m6A workflow, please follow the instructions provided in the original [GLORI-tools](https://github.com/liucongcas/GLORI-tools).

To prepare TSS annotations for the m6Am workflow, use:

```bash
# ${DUOdir} is the directory containing GLORI-DUO-tools scripts
# ${gtf2} is the GTF annotation file

python ${DUOdir}/Pre/gtf2tbl.py -i ${gtf2} -o ${gtf2}.anno

# Generate regions from 100 bp upstream to transcript ends
# (used for annotation of identified TSS sites)

cat ${gtf2}.anno | awk 'BEGIN {FS="\t";OFS="\t"} NR>1 {if ($3=="+") {print $2,$4-100,$5,$13,$1,$3,$4+1,$14,$15} else {print $2,$4,$5+100,$13,$1,$3,$5,$14,$15}}' > ${gtf2}.TSS.u100dInf.bed
```

---

# One-command Site Calling

Examples of reference and annotation files:

```bash
genome="${anno_dir}/hg38_chr_only.fa.AG_conversion.fa"
genome2="${anno_dir}/hg38_chr_only.fa"
rvsgenome="${anno_dir}/hg38_chr_only.rvsCom.fa"
TfGenome="${anno_dir}/GCF_000001405.39_GRCh38.p13_rna2.fa.AG_conversion.fa"

anno="${anno_dir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2"
baseanno="${anno_dir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base"

tssanno="${gtf2}.TSS.u100dInf.bed"

gtf="${anno_dir}/GCF_000001405.39_GRCh38.p13_genomic.AG_converted.gtf"
gtf2="${anno_dir}/GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens"

DUOdir=${path_to_GLORI-DUO-tools}
```

To identify m6A(m) sites from FASTQ files, use:

```bash
# for m6Am

python ${DUOdir}/DUO.py --raw_fq ${rawfq} -o ${outdir} \
  --mode m6Am --prx ${prx} \
  -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} \
  -a ${anno} -ta ${tssanno} -ba ${baseanno} --gtf ${gtf} --gtf2 ${gtf2} \
  >${log_out} 2>${log_err}

# Output files:
# ${outdir}/${prx}_TSS_raw.bed.annotated.rmdup.gz
# ${outdir}/${prx}_m6Am_sites_raw.tsv.gz
```

```bash
# for m6A

python ${DUOdir}/DUO.py --raw_fq ${rawfq} -o ${outdir} \
  --mode m6A --prx ${prx} \
  -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} \
  -a ${anno} -ta ${tssanno} -ba ${baseanno} --gtf ${gtf} --gtf2 ${gtf2} \
  >${log_out} 2>${log_err}

# Output file:
# ${outdir}/${prx}.totalm6A.FDR.csv.gz
```

The output files contain unfiltered m6A(m) sites, which can be further filtered and merged across multiple samples using:

```bash
# for m6Am

python ${DUOdir}/Call_m6Am/merge_multisam_m6Am.py \
      -i ${outdir}/*_TSS_raw.bed.annotated.rmdup.gz \
      -o ${merged_m6Am} \
      -c 15 -C 5 --tpm 0.5 --prop 0.05 \
      2>${log}
```

```bash
# for m6A

python ${DUOdir}/Call_m6A/merge_multisam_m6A.py \
      -i ${outdir}/${prx}.totalm6A.FDR.csv.gz \
      -o ${merged_m6A} \
      -c 15 -C 5 -s 0.8 -r 0.1 -adp 0.05 \
      2>${log}
```

---

# Maintainers and Contributing

GLORI-DUO-tools is developed and maintained by Zhe Zhou (zzhou24@pku.edu.cn).

Contributions, issue reports, and pull requests are welcome.

GitHub repository:

https://github.com/ZedekiahZhou/GLORI-DUO-tools

---

# License

Released under the MIT License.

---
