configfile: "config.yaml"

from snakemake.utils import makedirs

#################################
# author: Carolina Pita Barros  #
# carolina.pitabarros@wur.nl    #
# date: October 2021            #
#################################

if "OUTDIR" in config:
    print("\nSaving to " + config["OUTDIR"] + "\n")
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

pipeline = "nanopore-assembly" # replace with your pipeline's name


include: "rules/create_file_log.smk"

LONGREADS = config["LONGREADS"]
SHORTREADS = config["SHORTREADS"]
PREFIX = config["PREFIX"]
GENOME_SIZE = config["GENOME_SIZE"]

if "OUTDIR" in config:
    # print("\nSaving to " + config["OUTDIR"] + "\n")
    workdir: config["OUTDIR"]

makedirs("logs_slurm")


rule all:
    input:
        files_log,
        # expand("{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa", prefix = PREFIX),
        # expand("2_assembly/{prefix}/assembly.fasta", prefix=PREFIX)
        # expand("polish_{prefix}.done", prefix=PREFIX),
        # expand("busco_{prefix}.done", prefix=PREFIX)
        expand("results/longshot_{prefix}.vcf", prefix=PREFIX),
        expand("results/assembly_stats_{prefix}.txt",prefix=PREFIX),
        expand("busco_{prefix}.done", prefix=PREFIX),
        expand("variant_calling/{prefix}.vcf.stats",prefix=PREFIX)

rule trimming_adaptors:
    input:
        os.path.join(workflow.basedir,LONGREADS)
    output:
        "1_trimming/{prefix}.trimmed.fastq"
    message:
        'Rule {rule} processing'
    shell:
        'porechop -i {input} -o {output} -t 8 --discard_middle'

rule assemble_flye:
    input:
        rules.trimming_adaptors.output
    output:
        assembly = "2_assembly/{prefix}/assembly.fasta",
        # multiext("2_assembly/{prefix}/assembly_graph", ".gfa", ".gv"),
        # info ="2_assembly/{prefix}/assembly_info.txt"
    message:
        'Rule {rule} processing'
    params:
        # genome_size = GENOME_SIZE,
        outdir = "2_assembly/{prefix}"
    shell:
        "flye --nano-raw {input} --out-dir {params.outdir} --threads 16"

rule one_line_fasta:
    input:
        rules.assemble_flye.output.assembly
    output:
        "{prefix}_oneline.fa"
    message:
        'Rule {rule} processing'
    log:
        err = "logs_slurm/one_line_fasta_{prefix}.err",
    shell:
        'seqtk seq -l0 {input} > {output} 2> {log.err}'

LONGREADS_PATH = os.path.join(workflow.basedir,LONGREADS)
rule scaffolding_long_reads:
    input:
        draft = rules.one_line_fasta.output,
        reads =  LONGREADS_PATH
    output:
        "{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa"
    message:
        'Rule {rule} processing'
    params:
        size=GENOME_SIZE,
        draft = os.path.splitext(rules.one_line_fasta.output[0])[0],
        reads = os.path.splitext(os.path.splitext(os.path.basename(LONGREADS))[0])[0]
    shell:
        """
        longstitch ntLink-arks draft={params.draft} reads={params.reads} G={params.size}
        """


rule busco:
    input:
        ref = rules.scaffolding_long_reads.output
    output:
        "{prefix}_vertebrae/short_summary.specific.aves_odb10.oldturkey_aves.txt"
    message:
        "Rule {rule} processing"
    params:
        outdir = "{prefix}_vertebrae"
    shell:
        "busco -m genome -i {input.ref} -c 12 -o {params.outdir} -l vertebrata_odb10"

polca_ext = [".alignSorted.bam",".alignSorted.bam.bai", ".fa.fai", ".batches", ".names", ".report" , ".bwa.bwt",".bwa.pac",".bwa.ann",".bwa.amb",".bwa.sa"]
polca_ext_temp = [".sort.success", ".vc.success", ".report.success", ".unSorted.sam", ".index.success", ".fix.success", "map.success"]

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
    shell:
        """
        # module load bwa
        polca.sh -a {input.assembly} -r '{input.reads}' -t 16 -m 2G
        """
        
rule index_vcf:
    input:
        rules.polish_polca.output.vcf
    output:
        vcf = "variant_calling/{prefix}.vcf",
        idx = "variant_calling/{prefix}.vcf.tbi"
    message:
        "Rule {rule} processing"
    shell:
        """
module load bcftools
mv {input} {output.vcf}
tabix -p vcf {output.vcf}
        """

rule bcftools_stats:
    input:
        rules.index_vcf.output.vcf
    output:
        "variant_calling/{prefix}.vcf.stats"
    message:
        'Rule {rule} processing'
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
    shell:
        """
        python {params.script} {input} > {output}
        """

rule minimap2:
    input:
        longreads = LONGREADS_PATH,
        assembly = rules.polish_polca.output.assembly
    output:
        'mapped/{prefix}_longreads.mapped.bam'
    message:
        'Rule {rule} processing'
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
    shell:
        """
module load samtools
samtools sort -@ 16 -T sort_temp -o {output} {input}
samtools index {output}
        """

rule var_calling_longshot:
    input:
        bam = rules.sort_index_longreads.output,
        assembly = rules.polish_polca.output.assembly
    output:
        "results/longshot_{prefix}.vcf"
    message:
        'Rule {rule} processing'
    shell:
        """
longshot --bam {input.bam} --ref {input.assembly} --out {output}
        """


# rule run_busco:
#     input:
#         rules.scaffolding_long_reads.output
#     output:
#         directory("{prefix}_vertebrae"),
#         touch("busco_{prefix}.done")
#     log:
#         "logs_slurm/busco_{prefix}.log"
#     threads: 8
#     params:
#         mode="genome",
#         lineage="vertebrata_odb10",
#         downloads_path="resources/busco_downloads",
#         # # optional parameters
#         # extra=""
#     wrapper:
#         "0.78.0/bio/busco"