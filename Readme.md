ToDo:
1. reads5p pileup multiprocessing
2. m6A样本的所有点信息？怎么和merge结合
3. m6Am用于merge之前的文件名需要更改
4. 若A.bed_sorted中read name与*.sam不匹配，awk脚本会卡死！！！
5. 支持input为多个分散的fastq，可以省略合并步骤


All changes compared to the original GLORI-tools:

1. use awk to accelerate AG conversion
2. add support for bowtie2 mapping
3. STAR maaping: add "--readNameSeparator ' '" to support BGI read names (which contain '/')
4. Fix the FDR calculation steps
5. split the main steps, allowing starting from intermediate files