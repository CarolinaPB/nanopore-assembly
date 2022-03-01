configfile: "config.yaml"

from snakemake.utils import makedirs, linecount

#################################
# author: Carolina Pita Barros  #
# carolina.pitabarros@wur.nl    #
# date: October 2021            #
#################################

if "OUTDIR" in config:
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

pipeline = "nanopore-assembly"


include: "rules/create_file_log.smk"

LONGREADS = config["LONGREADS"]
SHORTREADS = config["SHORTREADS"]
PREFIX = config["PREFIX"]
GENOME_SIZE = config["GENOME_SIZE"]
BUSCO_LINEAGE = config["BUSCO_LINEAGE"]

# If COMPARISON_GENOME is set, compute whole genome alignment between polished assembly and comparison genome
if "COMPARISON_GENOME" in config:
    comp_genome_results = expand("results/genome_alignment/{prefix}_{species}.png", prefix=PREFIX, species = config["COMPARISON_GENOME"].keys())
    MIN_ALIGNMENT_LENGTH = config["MIN_ALIGNMENT_LENGTH"]
    MIN_QUERY_LENGTH = config["MIN_QUERY_LENGTH"]
else:
    comp_genome_results = []

localrules: longreads_softlink, create_file_log

rule all:
    input:
        files_log,
        expand("results/assembly_stats_{prefix}.txt",prefix=PREFIX),
        expand("results/variant_calling/{prefix}_shortreads.vcf.gz.stats",prefix=PREFIX),
        expand("results/variant_calling/{prefix}_longreads.vcf.gz.stats",prefix=PREFIX),
        expand("busco_{prefix}_{busco_cat}_polish_{lineage}/short_summary.specific.{lineage}.busco_{prefix}_{busco_cat}_polish_{lineage}.txt", prefix=PREFIX, lineage=BUSCO_LINEAGE, busco_cat = ["before", "after"]),
        comp_genome_results,


####################### ASSEMBLY #######################


LONGREADS_PATH = os.path.join(workflow.basedir,LONGREADS)

rule longreads_softlink:
    input:
        LONGREADS_PATH
    output:
        os.path.basename(LONGREADS)
    message:
        'Rule {rule} processing'
    shell:
        'ln -s {input} {output}'

rule trimming_adaptors:
    input:
        os.path.join(workflow.basedir,LONGREADS)
    output:
        "trimming/{prefix}.trimmed.fastq"
    message:
        'Rule {rule} processing'
    group:
        'assembly'
    shell:
        'porechop -i {input} -o {output} -t 8 --discard_middle'

rule assemble_flye:
    input:
        rules.trimming_adaptors.output
    output:
        assembly = "assembly/{prefix}/assembly.fasta"
    message:
        'Rule {rule} processing'
    params:
        # genome_size = GENOME_SIZE,
        outdir = "assembly/{prefix}"
    group:
        'assembly'
    shell:
        "flye --nano-raw {input} --out-dir {params.outdir} --threads 16"


####################### SCAFFOLDING #######################


rule one_line_fasta:
    input:
        rules.assemble_flye.output.assembly 
    output:
        "{prefix}_oneline.fa"
    message:
        'Rule {rule} processing'
    log:
        err = "logs_slurm/one_line_fasta_{prefix}.err",
    group:
        'scaffolding'
    shell:
        "seqtk seq -l0 {input} > {output} 2> {log.err}"

rule scaffolding_long_reads:
    input:
        draft = rules.one_line_fasta.output,
        reads =  rules.longreads_softlink.output
    output:
        "{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa"
    message:
        'Rule {rule} processing'
    params:
        size=GENOME_SIZE,
        draft = os.path.splitext(rules.one_line_fasta.output[0])[0],
        reads = os.path.splitext(os.path.splitext(os.path.basename(LONGREADS))[0])[0]
    group:
        'scaffolding'
    shell:
        """
        longstitch ntLink-arks draft={params.draft} reads={params.reads} G={params.size}
        """


####################### POLISHING #######################


polca_ext = [".alignSorted.bam",".alignSorted.bam.bai", ".fai", ".batches", ".names", ".report" , ".bwa.bwt",".bwa.pac",".bwa.ann",".bwa.amb",".bwa.sa"]
polca_ext_temp = [".sort.success", ".vc.success", ".report.success", ".unSorted.sam", ".index.success", ".fix.success", ".map.success"]

rule polish_polca:
    input:
        assembly = rules.scaffolding_long_reads.output,
        reads = SHORTREADS
    output:
        not_temp = expand("{{prefix}}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa{ext}", ext = polca_ext),
        to_remove =temp(expand("{{prefix}}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa{ext}", ext = polca_ext_temp)),
        assembly = "{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.PolcaCorrected.fa",
        vcf = "{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.vcf"
    message:
        'Rule {rule} processing'
    group:
        'polishing'
    shell:
        """
        # module load bwa
        polca.sh -a {input.assembly} -r '{input.reads}' -t 16 -m 2G
        """

rule index_polished_assembly:
    input:
        rules.polish_polca.output.assembly
    output:
        "{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.PolcaCorrected.fa.fai"
    message:
        'Rule {rule} processing'
    group:
        'polishing'
    shell:
        """
module load samtools
samtools faidx {input}
        """

####################### ASSEMBLY ASSESSMENT #######################


rule get_assembly_stats:
    input:
        rules.polish_polca.output.assembly,
    output:
        "results/assembly_stats_{prefix}.txt"
    message:
        'Rule {rule} processing'
    params:
        script = os.path.join(workflow.basedir, "scripts/get_assembly_stats.py"),
    group:
        'polishing'
    shell:
        """
        python {params.script} {input} > {output}
        """

def busco_input(wildcards):
    '''
    Function to give the correct input for busco before and after polishing
    '''
    if wildcards.busco_cat == "before":
        return(f"{wildcards.prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa")
    elif wildcards.busco_cat == "after":
        return(f"{wildcards.prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.PolcaCorrected.fa")

rule busco:
    input:
        busco_input
    output:
        "busco_{prefix}_{busco_cat}_polish_{lineage}/short_summary.specific.{lineage}.busco_{prefix}_{busco_cat}_polish_{lineage}.txt"
    message:
        'Rule {rule} processing'
    params:
        outdir = "busco_{prefix}_{busco_cat}_polish_{lineage}"
    shell:
        "busco -m genome -f -i {input} -c 12 -o {params.outdir} -l {wildcards.lineage}"



####################### VARIANT CALLING #######################


rule bwa_index:
    input: 
        rules.polish_polca.output.assembly
    output:
        multiext("{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.PolcaCorrected.fa", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
    group:
        "short_var_calling"
    conda:
        "envs/bwamem2.yaml"
    shell:
        "bwa-mem2 index {input}"

rule map_short_reads_bwa:
    input:
        assembly = rules.polish_polca.output.assembly,
        idx = rules.bwa_index.output,
        reads=SHORTREADS
    output:
        temp("mapped/{prefix}_shortreads.mapped.bam")
    resources: 
        cpus=16
    group:
        "short_var_calling"
    message:
        "Rule {rule} processing"
    conda:
        "envs/bwamem2.yaml"
    shell:
        """
        module load samtools
        bwa-mem2 mem -t {resources.cpus} {input.assembly} {input.reads} | samblaster -r | samtools view -b - > {output}
        """

rule map_longreads_minimap:
    input:
        longreads = LONGREADS_PATH,
        assembly = rules.polish_polca.output.assembly
    output:
        temp('mapped/{prefix}_longreads.mapped.bam')
    message:
        'Rule {rule} processing'
    shell:
        """
module load samtools 
minimap2 -t 16 --MD -a -x map-ont {input.assembly} {input.longreads}  | samtools view -S -b - > {output}
        """

rule sort_index_reads:
    input:
        'mapped/{prefix}_{read_type}.mapped.bam'
    output:
        bam = 'mapped/{prefix}_{read_type}.mapped.sorted.bam',
        bai = 'mapped/{prefix}_{read_type}.mapped.sorted.bam.bai'
    message:
        'Rule {rule} processing'
    shell:
        """
module load samtools
samtools sort -@ 16 -T sort_temp -o {output.bam} {input}
samtools index {output.bam}
        """

rule var_calling_freebayes:
    input:
        ref=rules.polish_polca.output.assembly,
        bam='mapped/{prefix}_shortreads.mapped.sorted.bam',
        indexes='mapped/{prefix}_shortreads.mapped.sorted.bam.bai',
        idx = rules.index_polished_assembly.output
    output:
        vcf = "results/variant_calling/{prefix}_shortreads.vcf.gz",
        idx = "results/variant_calling/{prefix}_shortreads.vcf.gz.tbi",
    params:
        chunksize=100000, # reference genome chunk size for parallelization (default: 100000)
        scripts_dir = os.path.join(workflow.basedir, "scripts")
    shell:
        """
module load freebayes bcftools vcflib python/2.7.15 samtools

{params.scripts_dir}/freebayes-parallel.sh <({params.scripts_dir}/fasta_generate_regions.py {input.ref}.fai {params.chunksize}) 2 \
-f {input.ref} \
--use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 \
{input.bam} | vcffilter -f 'QUAL > 20' | bgzip -c > {output.vcf}
tabix -p vcf {output.vcf}
        """


rule var_calling_longshot:
    input:
        bam = 'mapped/{prefix}_longreads.mapped.sorted.bam',
        assembly = rules.polish_polca.output.assembly,
        idx = rules.index_polished_assembly.output
    output:
        temp("tmp_var_calling/{prefix}_longreads.vcf")
    message:
        'Rule {rule} processing'
    group:
        'var_calling'
    shell:
        """
longshot --bam {input.bam} --ref {input.assembly} --out {output}
        """

rule filter_index_vcf:
    input:
        rules.var_calling_longshot.output
    output:
        vcf = "results/variant_calling/{prefix}_longreads.vcf.gz",
        idx = "results/variant_calling/{prefix}_longreads.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    shell:
        """
module load samtools vcflib/gcc/64/0.00.2019.07.10
vcffilter -f 'QUAL > 20' {input} | bgzip -c > {output.vcf}
tabix -p vcf {output.vcf}
        """

rule bcftools_stats:
    input:
        "results/variant_calling/{prefix}_{read_type}.vcf.gz"
    output:
        "results/variant_calling/{prefix}_{read_type}.vcf.gz.stats"
    message:
        'Rule {rule} processing'
    shell:
        """
module load bcftools
bcftools stats {input} > {output}
        """


####################### WHOLE GENOME ALIGNMENT #######################



def get_ref_path(wildcards):
    '''
    Get genome path for comparison species
    '''
    return(config["COMPARISON_GENOME"][wildcards.species])


rule align_genomes:
    input:
        assembly = rules.polish_polca.output.assembly,
        comparison = get_ref_path
    output:
        "results/genome_alignment/{prefix}_vs_{species}.paf"
    message:
        'Rule {rule} processing'
    group:
        'genome_alignment'
    shell:
        """
minimap2 -t 12 -cx asm5 {input.comparison} {input.assembly} > {output}
        """

rule plot_aligned_genomes:
    input:
        rules.align_genomes.output
    output:
        "results/genome_alignment/{prefix}_{species}.png"
    message:
        'Rule {rule} processing'
    params:
        script = os.path.join(workflow.basedir, "scripts/pafCoordsDotPlotly.R"),
        min_alignment_length = MIN_ALIGNMENT_LENGTH,
        min_query_length = MIN_QUERY_LENGTH, 
        outdir = "results/genome_alignment/"
    group:
        'genome_alignment'
    shell:
        """
module load R
Rscript {params.script} -i {input} -o {wildcards.prefix}_{wildcards.species} -s -t -x -m {params.min_alignment_length} -q {params.min_query_length} -l
mv {wildcards.prefix}_{wildcards.species}.png {params.outdir}
        """


####################### CLEAN UP #######################


onsuccess:
    print("Workflow finished, no error")
    shell("mv *oneline*PolcaCorrected* results")
    shell("mv *log logs_slurm")
    shell("mkdir -p other_files")
    shell("mv *oneline* other_files")
