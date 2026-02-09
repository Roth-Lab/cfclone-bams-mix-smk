from snakemake.utils import min_version, validate

min_version("9.12")

from utils import ConfigManager

config = ConfigManager(config)

rule all:
    input:
        config.pipeline_files

onsuccess:
    config.email_notification(
        on="success",
        workflow="cfclone-bams-mix-tf-smk",
        configfile=workflow.configfiles[0],
        imgs=[config.copied_config],
    )

onerror:
    config.email_notification(
        on="error", 
        workflow="cfclone-bams-mix-tf-smk",
        configfile=workflow.configfiles[0],
    )

pathvars:
    out_dir=str(config.replicate_out_dir),
    pipeline_dir=str(config.replicate_pipeline_dir),


ruleorder: down_sample_bam_file > merge_down_sampled_bam_files > write_down_sampled_bam_summary_file > write_mixed_bam_summary_file

rule build_config_file:
    input:
        workflow.configfiles[0]
    output:
        config.copied_config
    shell:
        "cp {input} {output}"


# DOWN SAMPLE BAM FILES 

rule down_sample_bam_file:
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
        config.compute_down_sample_proportion
    shell:
        """
        samtools view -b -s {params} {input.bam} -o {output.bam} 1> {log} 2>&1
        samtools index {output.bam} 1>> {log} 2>&1
        """
        # samtools view --write-index -b -s {params.fraction} {input.bam} -o {output.bam}##idx##{output.bai} 1> {log} 2>&1 # returns error "Random alignment retrieval only works for index ...""

rule write_down_sampled_bam_summary_file:
    input:
        bam=config.down_sampled_bam_file_template,
        bai=config.down_sampled_bai_file_template,
    output:
        config.down_sampled_total_reads_template
    conda:
        "envs/python.yaml"
    threads: 2
    log:
        config.get_log_file(config.down_sampled_total_reads_template)
    params:
        r=config.read_length,
        g=config.genome_length,
        p=config.patient,
        s=config.get_sample, 
        c=lambda wildcards: config.coverages[int(wildcards.coverage_id)],
        i=lambda wildcards: config.proportions[int(wildcards.proportion_id)],
    shell:
        "(python scripts/write_down_sample_summary_file.py "
        "-i {input.bam} "
        "-o {output} "
        "--patient {params.p} "
        "--sample {params.s} "
        "--coverage {params.c} "
        "--proportion {params.i} "
        "--bam-id {wildcards.bam_id} "
        "--read-length {params.r} "
        "--genome-length {params.g} ) >{log} 2>&1"


rule merge_down_sampled_bam_files:
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
        r=config.read_length,
        g=config.genome_length,
        p=config.patient,
        s=config.get_sample, 
        c=lambda wildcards: config.coverages[int(wildcards.coverage_id)],
        i=lambda wildcards: config.proportions[int(wildcards.proportion_id)],
    shell:
        "(python scripts/write_down_sample_summary_file.py "
        "-i {input.bam} "
        "-o {output} "
        "--patient {params.p} "
        "--sample {params.s} "
        "--coverage {params.c} "
        "--proportion {params.i} "
        "--bam-id 'mixed' "
        "--read-length {params.r} "
        "--genome-length {params.g} ) >{log} 2>&1"


rule merge_down_sampled_summary_files:
    input:
        config.gather_files(config.down_sampled_total_reads_template)
    output:
        config.down_sampled_summary_file
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.down_sampled_summary_file)
    shell:
        "(python scripts/merge_tables.py -i {input} -o {output} ) >{log} 2>&1"


rule merge_mixed_bams_summary_files:
    input:
        config.gather_mixed_bams_summary_files
    output:
        config.mixed_bams_summary_file
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.mixed_bams_summary_file)
    shell:
        "(python scripts/merge_tables.py -i {input} -o {output} ) >{log} 2>&1"


# # BUILD RDR AND BAF TSV FILES 

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
        "../cfclone-smk/Snakefile"
    config:
        config.cfclone_config


use rule * from cfclone as cfclone_*


use rule run_cfclone from cfclone as cfclone_run_cfclone with:
    input:
        c=config.cfclone_clone_cn_template,
        i=config.cfclone_ctdna_template


rule build_replicate_summary:
    input:
        e=config.replicate_evidence_template,
        t=config.replicate_tumour_content_template,
    output:
        config.replicate_summary_template
    params:
        p=config.patient,
        s=config.get_sample, 
        c=lambda wildcards: config.coverages[int(wildcards.coverage_id)],
        i=lambda wildcards: config.proportions[int(wildcards.proportion_id)],
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.replicate_summary_template)
    shell:
        "(python scripts/write_summary_file.py "
        "-e {input.e} "
        "-t {input.t} "
        "-o {output} "
        "--patient {params.p} "
        "--sample {params.s} "
        "--coverage {params.c} "
        "--proportion {params.i}) >{log} 2>&1"


rule merge_summaries:
    input:
        config.gather_cfclone_summary_files
    output:
        config.summary_file
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.summary_file),
    shell:
        "(python scripts/merge_tables.py -i {input} -o {output}) >{log} 2>&1"


rule plot_summaries:
    input:
        config.summary_file,
    output:
        config.plot_summary_file,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.plot_summary_file)
    shell:
        "(python scripts/plot_tf.py -i {input} -o {output}) >{log} 2>&1"