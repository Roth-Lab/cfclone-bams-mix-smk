


### Requirements before running workflow
1. Needs directory that contains bams files.
2. Need .tsv file that contains read counts for each bam file in same directory as bam files.
3. See `resources/counts_bams.sh`

### Requirements to run workflow




### Notes

- All files and directories below can be templated using {patient_id} to allow for config file reuse e.g. snp_file: /foo/bar/data/{patient_id}.bcf

- control bam file ideally should be ctDNA sequencing from a copy number flat genome. Not the matached normal for the patient.

- Path to file with SNPs for the patient either BCF or VCF format
- Note: SNPs should be called from the matched normal
- Note: SNPs should be phased