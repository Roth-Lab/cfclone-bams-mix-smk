from snakemake.utils import min_version, validate

min_version("9.12")

validate(config, "schemas/config.schema.yaml")

from utils import ConfigManager

config = ConfigManager(config)


pathvars:
    out_dir=str(config.cfclone_out_dir),
    pipeline_dir=str(config.cfclone_pipeline_dir),

ruleorder: downsample_bam_file > merge_bam_files > write_mixed_bam_summary_file

rule all:
    input:
        config.pipeline_files

rule build_config_file:
    input:
        workflow.configfiles[0]
    output:
        config.copied_config
    shell:
        "cp {input} {output}"

# DOWN SAMPLE BAM FILES 

rule downsample_bam_file:
    input:
        bam=config.get_bam_file,
        bai=config.get_bai_file,
    output:
        bam=temp(config.down_sampled_bam_file_template),
        idx=temp(config.down_sampled_bai_file_template),
    conda:
        "envs/samtools.yaml"
    threads: 
        2
    log:
        config.get_log_file(config.down_sampled_bam_file_template)
    benchmark:
        config.get_benchmark_file(config.down_sampled_bam_file_template)
    params:
        prop=config.compute_down_sample_proportion,
        seed=lambda wildcards: int(wildcards.data_seed_id)
    shell:
        """
        samtools view -b -f 0x2 --subsample-seed {params.seed} --subsample {params.prop} {input.bam} -o {output.bam} 1> {log} 2>&1
        samtools index {output.bam} 1>> {log} 2>&1
        """
        # samtools view -b -s {params} {input.bam} -o {output.bam} 1> {log} 2>&1
        # samtools view --write-index -b -s {params.fraction} {input.bam} -o {output.bam}##idx##{output.bai} 1> {log} 2>&1 # returns error "Random alignment retrieval only works for index ...""


rule merge_bam_files:
    input:
        config.get_bams_to_mix
    output:
        bam=temp(config.mixed_bam_file_template),
        bai=temp(config.mixed_bai_file_template)
    log:
        config.get_log_file(config.mixed_bam_file_template)
    threads: 8
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools merge {output.bam} {input} 1> {log} 2>&1
        samtools index {output.bam} 1>> {log} 2>&1
        """


rule write_mixed_bam_summary_file:
    input:
        bam=config.mixed_bam_file_template,
        bai=config.mixed_bai_file_template,
    output:
        config.mixed_bam_total_reads_template
    conda:
        "envs/python.yaml"
    threads: 2
    log:
        config.get_log_file(config.mixed_bam_total_reads_template)
    params:
        p=config.patient,
        c=lambda wildcards: config.coverages[int(wildcards.coverage_id)],
        s=config.get_mixed_samples, 
        i=lambda wildcards: config.proportions[int(wildcards.proportion_id)],
        d=lambda wildcards: int(wildcards.data_seed_id),
        r=config.read_length,
        g=config.genome_length
    shell:
        "(python scripts/write_bam_summary_file.py "
        "-i {input.bam} "
        "-o {output} "
        "--patient {params.p} "
        "--coverage-id {wildcards.coverage_id} "
        "--mixture-id {wildcards.mixture_id} "
        "--proportion-id {wildcards.proportion_id} "
        "--data-seed {params.d} "
        "--bam-id 'mixed' "
        "--coverage {params.c} "
        "--mixture {params.s} "
        "--proportion {params.i} "
        "--read-length {params.r} "
        "--genome-length {params.g} ) >{log} 2>&1"

rule merge_mixed_bams_summary_files:
    input:
        config.gather_files(config.mixed_bam_total_reads_template)
    output:
        config.mixed_bams_summary_file
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.mixed_bams_summary_file)
    shell:
        "(python scripts/merge_tables.py -i {input} -o {output} ) >{log} 2>&1"


# BUILD RDR AND BAF TSV FILES 

rule build_map_wig:
    input:
        config.map_file
    output:
        temp(config.map_template)
    params:
        c=",".join(config.chromosomes),
        w=config.parse_region_size(config.bin_size)
    conda:
        "envs/hmmcopy-utils.yaml"
    log:
        config.get_log_file(config.map_template)
    resources:
        mem="8G"
    shell:
        "mapCounter -c {params.c} -w {params.w} -s {input} "
        "> {output} "
        "2> {log}"

rule build_gc_wig:
    input:
        config.ref_genome_file
    output:
        temp(config.gc_template)
    params:
        c=",".join(config.chromosomes),
        w=config.parse_region_size(config.bin_size)
    conda:
        "envs/hmmcopy-utils.yaml"
    log:
        config.get_log_file(config.gc_template)
    resources:
        mem="8G"
    shell:
        "gcCounter -c {params.c} -w {params.w} -s {input} "
        "> {output} "
        "2>{log}"

rule build_reads_chrom:
    input:
        bam=config.mixed_bam_file_template,
        bai=config.mixed_bai_file_template
    output:
        temp(config.reads_chrom_template)
    params:
        q=config.min_mqual,
        s=config.parse_region_size(config.bin_size),
    conda:
        "envs/python.yaml"
    wildcard_constraints:
        chrom='[^/]+'
    log:
        config.get_log_file(config.reads_chrom_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/build_read_counts.py "
        "-b {input.bam} "
        "-o {output} "
        "-c {wildcards.chrom} "
        "-q {params.q} "
        "-s {params.s}) >{log} 2>&1"


rule build_normal_reads_chrom:
    input:
        bam=config.get_control_bam_file,
        bai=config.get_control_bai_file
    output:
        temp(config.normal_reads_chrom_template)
    params:
        q=config.min_mqual,
        s=config.parse_region_size(config.bin_size),
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.normal_reads_chrom_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/build_read_counts.py "
        "-b {input.bam} "
        "-o {output} "
        "-c {wildcards.chrom} "
        "-q {params.q} "
        "-s {params.s}) >{log} 2>&1"
        

rule build_reads:
    input:
        lambda wildcards: config.gather_files_by_chrom(config.reads_chrom_template, wildcards)
    output:
        temp(config.reads_template)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.reads_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/merge_tables.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"


rule build_normal_reads:
    input:
        lambda wildcards: config.gather_files_by_chrom(config.normal_reads_chrom_template, wildcards)
    output:
        temp(config.normal_reads_template)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.normal_reads_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/merge_tables.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"


rule build_reads_wig:
    input:
        config.reads_template
    output:
        config.reads_wig_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.reads_wig_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/build_reads_wig.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"

rule build_rdr:
    input:
        i=config.reads_template,
        b=config.black_list_file,
        c=config.centromere_file,
        g=config.gc_template,
        m=config.map_template,
        n=config.normal_reads_template
    output:
        temp(config.rdr_template)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.rdr_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/build_corrected_rdr.py "
        "-i {input.i} "
        "-b {input.b} "
        "-c {input.c} "
        "-g {input.g} "
        "-m {input.m} "
        "-n {input.n} "
        "-o {output}) >{log} 2>&1"


rule build_allele_counts_chrom:
    input:
        bam=config.mixed_bam_file_template,
        bai=config.mixed_bam_file_template,
        s=config.get_snp_file
    output:
        temp(config.allele_counts_chrom_template)
    params:
        q=config.min_bqual,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.allele_counts_chrom_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/build_allele_counts.py "
        "-b {input.bam} "
        "-s {input.s} "
        "-o {output} "
        "-c {wildcards.chrom} "
        "-q {params.q}) >{log} 2>&1"


rule build_allele_counts:
    input:
        lambda wildcards: config.gather_files_by_chrom(config.allele_counts_chrom_template, wildcards)
    output:
        temp(config.allele_counts_template)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.allele_counts_template)
    shell:
        "(python scripts/merge_tables.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"


rule build_hap_bin_counts:
    input:
        bam=config.mixed_bam_file_template,
        bai=config.mixed_bai_file_template,
        i=config.allele_counts_template
    output:
        temp(config.baf_template)
    params:
        c=" ".join(config.chromosomes),
        s=config.parse_region_size(config.bin_size),
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.baf_template)
    shell:
        "(python scripts/build_hap_bin_counts.py "
        "-b {input.bam} "
        "-i {input.i} "
        "-o {output} "
        "-c {params.c} "
        "-s {params.s} ) >{log} 2>&1"


rule build_combined_results:
    input:
        b=config.baf_template,
        r=config.rdr_template
    output:
        config.combined_results_template
    params:
        c=" ".join(config.chromosomes),
        s=config.parse_region_size(config.bin_size),
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.combined_results_template)
    shell:
        "(python scripts/build_combined_results.py "
        "-b {input.b} "
        "-r {input.r} "
        "-o {output} ) >{log} 2>&1"


rule plot_baf:
    input:
        config.combined_results_template
    output:
        config.baf_plot_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.baf_plot_template)
    shell:
        "(python scripts/plot_baf.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"


rule plot_rdr:
    input:
        i=config.combined_results_template
    output:
        config.rdr_plot_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.rdr_plot_template)
    shell:
        "(python scripts/plot_rdr.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"


# RUN CFCLONE 


rule build_cfclone_clone_file:
    input:
        c=config.get_clone_filter_file,
        i=config.get_hapclone_results_file,
    output:
        config.cfclone_clone_cn_template,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.cfclone_clone_cn_template),
    shell:
        "(python scripts/build_clone_cn_file.py -c {input.c} -i {input.i} -o {output}) >{log} 2>&1"


rule build_cfclone_ctdna_file:
    input:
        config.combined_results_template,
    output:
        config.cfclone_ctdna_template,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.cfclone_ctdna_template),
    shell:
        "(python scripts/build_ctdna_file.py -i {input} -o {output}) >{log} 2>&1"


module cfclone:
    snakefile:
        # "../cfclone-smk/Snakefile"
        "../cfclone-smk/Snakefile"
    config:
        config.cfclone_config


use rule * from cfclone as cfclone_*


use rule run_cfclone from cfclone as cfclone_run_cfclone with:
    input:
        c=config.cfclone_clone_cn_template,
        i=config.cfclone_ctdna_template,


rule build_summary_file:
    input:
        e=config.merged_evidence_file,
        t=config.merged_tumour_content_file,
    output:
        config.cfclone_summary_file
    params:
        p=config.patient,
        i=config.get_initial_bam_id,
        f=config.get_final_bam_id,
        c=config.get_coverage,
        t=config.get_proportion
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.cfclone_summary_file)
    shell:
        "(python scripts/write_summary_file.py "
        "--out-file {output} "
        "--evidence-file {input.e} "
        "--tumour-content-file {input.t} "
        "--patient {params.p} "
        "--initial-sample {params.i} "
        "--final-sample {params.f} "
        "--coverage {params.c} "
        "--proportion {params.t} "
        "--data-seed {wildcards.data_seed_id} ) >{log} 2>&1"


rule merge_summaries:
    input:
        config.gather_files(config.cfclone_summary_file)
    output:
        config.summary_file,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.summary_file),
    shell:
        "(python scripts/merge_tables.py -i {input} -o {output}) >{log} 2>&1"
