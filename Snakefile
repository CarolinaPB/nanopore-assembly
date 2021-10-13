configfile: "config.yaml"

from snakemake.utils import makedirs

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


rule all:
    input:
        files_log,
        expand("results/assembly_stats_{prefix}.txt",prefix=PREFIX),
        expand("variant_calling/{prefix}_shortreads.vcf.gz.stats",prefix=PREFIX),
        expand("variant_calling/{prefix}_longreads.vcf.gz.stats",prefix=PREFIX),
        expand("busco_{prefix}_scaffolded_polished/short_summary.specific.{lineage}.{prefix}_scaffolded_polished.txt", prefix=PREFIX, lineage=BUSCO_LINEAGE),
        expand("busco_{prefix}_scaffolded/short_summary.specific.{lineage}.{prefix}_scaffolded.txt", prefix=PREFIX, lineage=BUSCO_LINEAGE)



rule trimming_adaptors:
    input:
        os.path.join(workflow.basedir,LONGREADS)
    output:
        temp("1_trimming/{prefix}.trimmed.fastq")
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
        assembly = temp("2_assembly/{prefix}/assembly.fasta"),

    message:
        'Rule {rule} processing'
    params:
        # genome_size = GENOME_SIZE,
        outdir = "2_assembly/{prefix}"
    group:
        'assembly'
    shell:
        "flye --nano-raw {input} --out-dir {params.outdir} --threads 16"

rule one_line_fasta:
    input:
        rules.assemble_flye.output.assembly 
    output:
        temp("{prefix}_oneline.fa")
    message:
        'Rule {rule} processing'
    log:
        err = "logs_slurm/one_line_fasta_{prefix}.err",
    group:
        'scaffolding'
    shell:
        'seqtk seq -l0 {input} > {output} 2> {log.err}'

LONGREADS_PATH = os.path.join(workflow.basedir,LONGREADS)

rule scaffolding_long_reads:
    input:
        draft = rules.one_line_fasta.output,
        reads =  LONGREADS_PATH
    output:
        temp("{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa") #make temp
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

rule busco_before_polish:
    input:
        rules.scaffolding_long_reads.output
    output:
        "busco_{prefix}_scaffolded/short_summary.specific.{lineage}.{prefix}_scaffolded.txt"
    message:
        "Rule {rule} processing"
    params:
        outdir = "busco_{prefix}_scaffolded",
    shell:
        "busco -m genome -f -i {input} -c 12 -o {params.outdir} -l {wildcards.lineage}"

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
        vcf = temp("{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.vcf") #make temp
    message:
        'Rule {rule} processing'
    shell:
        """
        # module load bwa
        polca.sh -a {input.assembly} -r '{input.reads}' -t 16 -m 2G
        """

rule busco_after_polish:
    input:
        rules.polish_polca.output.assembly
    output:
        "busco_{prefix}_scaffolded_polished/short_summary.specific.{lineage}.{prefix}_scaffolded_polished.txt"
    message:
        "Rule {rule} processing"
    params:
        outdir = "busco_{prefix}_scaffolded_polished"
    shell:
        "busco -m genome -f -i {input} -c 12 -o {params.outdir} -l vertebrata_odb10"

rule index_vcf_shortreads:
    input:
        rules.polish_polca.output.vcf
    output:
        vcf = "variant_calling/{prefix}_shortreads.vcf.gz",
        idx = "variant_calling/{prefix}_shortreads.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    group:
        'postpolishing'
    shell:
        """
module load samtools
bgzip -c {input} > {output.vcf}
tabix -p vcf {output.vcf}
        """


rule bcftools_stats:
    input:
        "variant_calling/{prefix}{type}.vcf.gz"
    output:
        "variant_calling/{prefix}{type}.vcf.gz.stats"
    message:
        'Rule {rule} processing'
    group:
        'postpolishing'
    shell:
        """
module load bcftools
bcftools stats {input} > {output}
        """

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
        'postpolishing'
    shell:
        """
        python {params.script} {input} > {output}
        """

rule minimap2:
    input:
        longreads = LONGREADS_PATH,
        assembly = rules.polish_polca.output.assembly
    output:
        temp('mapped/{prefix}_longreads.mapped.bam')
    message:
        'Rule {rule} processing'
    group:
        'var_calling'
    shell:
        """
module load samtools 
minimap2 -t 16 --MD -a -x map-ont {input.assembly} {input.longreads}  | samtools view -S -b - > {output}
        """

rule sort_index_longreads:
    input:
        rules.minimap2.output
    output:
        'mapped/{prefix}_longreads.mapped.sorted.bam'
    message:
        'Rule {rule} processing'
    group:
        'var_calling'
    shell:
        """
module load samtools
samtools sort -@ 16 -T sort_temp -o {output} {input}
samtools index {output}
        """

rule index_polished_assembly:
    input:
        rules.polish_polca.output.assembly
    output:
        "{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.PolcaCorrected.fa.fai"
    message:
        'Rule {rule} processing'
    group:
        'postpolishing'
    shell:
        """
module load samtools
samtools faidx {input}
        """

rule var_calling_longshot:
    input:
        bam = rules.sort_index_longreads.output,
        assembly = rules.polish_polca.output.assembly,
        idx = rules.index_polished_assembly.output
    output:
        temp("longshot_{prefix}.vcf")
    message:
        'Rule {rule} processing'
    group:
        'var_calling'
    shell:
        """
longshot --bam {input.bam} --ref {input.assembly} --out {output}
        """

rule index_vcf_longshot:
    input:
        rules.var_calling_longshot.output
    output:
        vcf = "variant_calling/{prefix}_longreads.vcf.gz",
        idx = "variant_calling/{prefix}_longreads.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    group:
        'var_calling'
    shell:
        """
module load samtools
bgzip -c {input} > {output.vcf}
tabix -p vcf {output.vcf}
        """
