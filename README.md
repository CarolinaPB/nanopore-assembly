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
- [Minimap2](https://github.com/lh3/minimap2) - map long reads to reference. Genome alignment
- [Samtools](http://www.htslib.org/) - sort and index mapped reads and vcf files
- [Longshot](https://github.com/pjedge/longshot) - variant calling with nanopore reads
- [bcftools](https://samtools.github.io/bcftools/bcftools.html) - vcf statistics
- R - [pafCoordsDotPlotly](https://github.com/tpoorten/dotPlotly) - plot genome alignment 

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
BUSCO_LINEAGE:
  - <lineage>

# genome alignment parameters:
COMPARISON_GENOME: 
  <species>: /path/to/genome/fasta

# filter alignments less than cutoff X bp
MIN_ALIGNMENT_LENGTH: 10000
MIN_QUERY_LENGTH: 50000
```
- LONGREADS - name of file with long reads. This file should be in the working directory (where this config and the Snakefile are)
- SHORTREADS - paths to short reads fq.gz
- GENOME_SIZE - approximate genome size ```haploid genome size (bp)(e.g. '3e9' for human genome)``` from [longstitch](https://github.com/bcgsc/longstitch#full-help-page)
- PREFIX -  prefix for the created files
- OUTDIR - directory where snakemake will run and where the results will be written to  
  If you want the results to be written to this directory (not to a new directory), open config.yaml and comment out `OUTDIR: /path/to/outdir`
- BUSCO_LINEAGE - lineage used for busco. Can be one or more (one per line). To see available lineages run `busco --list-datasets`
- COMPARISON_GENOME - genome for whole genome comparison. Add your species name and the path to the fasta file. ex: `chicken: /path/to/chicken.fna.gz`. You can add several genomes, one on each line.   
  - If you don't want to run the genome alignment step, comment out 
```
COMPARISON_GENOME: 
  <species>: /path/to/genome/fasta
```
- MIN_ALIGNMENT_LENGTH and MIN_QUERY_LENGTH - parameters for plotting. If your plot is coming out blank or if there's an error with the plotting step, try lowering these thresholds. This happens because the alignments are not large enough.


If you have your long reads in several fastq files and need to create one file compressed file with all the reads:
1. In your pipeline directory create one file with all the reads
```
cat /path/to/fastq/directory/*.fastq > <name of file>.fq
```
2. Compress the file you just created:
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
The working directory will be messy with all the necessary files and results from the several pipeline steps.
The most important files are and directories are:  
- **<run_date>_files.txt** dated file with an overview of the files used to run the pipeline (for documentation purposes)
- **results** directory that contains
  - assembly_stats_\<prefix>.txt file with assembly statistics for the final assembly
  - **busco_{prefix}__before_polish_** and **busco_{prefix}_after_polish** directories - contain busco results before and after polishing respectively
    - short_summary.specific.{lineage}.{prefix}_before_polish.txt
    - short_summary.specific.{lineage}.{prefix}_after_polish.txt"
  - **variant_calling** directory with variant calling VCF files with long and short reads, as well as VCF stats
    - {prefix}_shortreads.vcf.gz
    - {prefix}_shortreads.vcf.gz.stats
    - {prefix}_longreads.vcf.gz
    - {prefix}_longreads.vcf.gz.stats
  - **3_mapped**
    - {prefix}_longreads.mapped.sorted.bam - long reads mapped to the new assembly

