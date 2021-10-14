# Assemble nanopore reads and do variant calling with short and long reads

## First follow the instructions here:
[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT
This is a pipeline that uses `Flye` to create a nanopore assembly. It also does variant calling with long and short reads.  
The pipeline starts by using `porechop` to trim the adaptors, then it uses `Flye` to create the assembly. After that, `ntLink-arks` from `Lonstitch` is used to scaffold the assembly using the nanopore reads. The scaffolded assembly is polished with `polca`. `Polca` also does variant calling with the short reads, while `longshot` does variant calling with the nanopore reads. To run `longshot`, first the long reads are aligned to the assembly with `minimap2`.  
In the end, in addition to your assembly and variant calling results, you'll also get assembly statistics and busco scores before and after the polishing.

#### Tools used:
- [Porechop](https://github.com/rrwick/Porechop) - trim adaptors
- [Flye](https://github.com/fenderglass/Flye) - assembly
- [Seqtk](https://github.com/lh3/seqtk) - convert fasta to one line fasta
- [LongStitch (ntLink-arks)](https://github.com/bcgsc/longstitch) - scaffolding with nanopore reads
- [BUSCO](https://busco.ezlab.org/) - assess assembly completeness
- [MaSuRCA (polca)](https://github.com/alekseyzimin/masurca) - polish assembly and do variant calling with short reads
- Python - get assembly stats
- [Minimap2](https://github.com/lh3/minimap2) - map long reads to reference
- [Samtools](http://www.htslib.org/) - sort and index mapped reads and vcf files
- [Longshot](https://github.com/pjedge/longshot) - variant calling with nanopore reads
- [bcftools](https://samtools.github.io/bcftools/bcftools.html) - vcf statistics


| ![DAG](https://github.com/CarolinaPB/nanopore-assembly/blob/master/workflow.png) |
|:--:|
|*Pipeline workflow* |


### Edit config.yaml with the paths to your files
```
LONGREADS: <nanopore_reads.fq.gz>
SHORTREADS:
  - /path/to/short/reads_1.fq.gz
  - /path/to/short/reads_2.fq.gz
GENOME_SIZE: <approximate genome size>
PREFIX: <prefix>
OUTDIR: /path/to/outdir
```
- LONGREADS - name of file with long reads. This file should be in the working directory (where this config and the Snakefile are)
- SHORTREADS - paths to short reads fq.gz
- GENOME_SIZE - approximate genome size ```haploid genome size (bp)(e.g. '3e9' for human genome)``` from [longstitch](https://github.com/bcgsc/longstitch#full-help-page)
- PREFIX -  prefix for the created files
- OUTDIR - directory where snakemake will run and where the results will be written to

If you have your long reads in in several fastq files and need to create one file compressed file with all the reads:
1. In your pipeline directory create one file with all the reads
```
cat /path/to/fastq/directory/*.fastq > <name of file>.fq
```
2. Compress the file you just greated:
```
gzip <name of file>.fq
```

### After installing the conda environment (first step of this guide) you'll need to edit the polca.sh file.
First go to the directory where miniconda3 is installed (usually your home directory). Go to `/<home>/miniconda/envs/<env_name>/bin` and open the file `polca.sh`. In my case the path looks like this: `/home/WUR/<username>/miniconda3/envs/<env_name>/bin/`. In your editor open `polca.sh` and replace this line:
```
$SAMTOOLS sort -m $MEM -@ $NUM_THREADS <(samtools view -uhS $BASM.unSorted.sam) $BASM.alignSorted 2>>samtools.err && \
```
With this:
```
$SAMTOOLS sort -m $MEM -@ $NUM_THREADS <(samtools view -uhS $BASM.unSorted.sam) -o $BASM.alignSorted.bam 2>>samtools.err && \
```

## RESULTS

