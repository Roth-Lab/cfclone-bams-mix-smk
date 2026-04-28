import pandas as pd

import numpy as np

from pathlib import Path

from snakemake.shell import shell


class ConfigManager:
    def __init__(self, config: dict):
        self.config = config

    # WILDCARDS

    @property
    def final_samples(self) -> list[str]:
        return self.config["final_samples"]

    @property
    def final_sample_ids(self) -> list[int]:
        return list(range(len(self.final_samples)))
    
    @property
    def mixtures(self) -> dict[int, dict[str, str]]:
        return self.config['mixtures']
    
    @property
    def mixture_ids(self) -> list[int]:
        return list(self.mixtures.keys())

    @property
    def increment(self) -> float:
        return self.config["increment"]

    @property
    def proportions(self) -> list[float]:
        return list(np.arange(0.0, 1.0 + self.increment, self.increment))

    @property
    def proportion_ids(self) -> list[int]:
        return list(range(len(self.proportions)))

    @property
    def coverages(self) -> list[float]:
        return self.config["coverages"]

    @property
    def coverage_ids(self) -> list[int]:
        return list(range(len(self.coverages)))
    
    @property
    def num_data_replicates(self) -> int:
        return self.config['num_data_replicates']
    
    @property
    def data_seed_ids(self) -> list[int]:
        return list(range(self.num_data_replicates))

    @property
    def bam_ids(self) -> list[str]:
        return ["initial", "final"]
    
    @property
    def num_model_replicates(self) -> int:
        return self.config['num_model_replicates']

    # DATA SETTINGS

    @property
    def patient(self) -> list[str]:
        return self.config["patient"]

    @property
    def initial_sample(self) -> str:
        return self.config["initial_sample"]

    @property
    def final_sample(self) -> str:
        return self.config["final_sample"]

    # INPUT COHORT DATA

    @property
    def cohort_tc_summary_file(self) -> str | None:
        return self.config.get('cohort_tc_summary_file', None)

    # INPUT PATIENT DATA

    @property
    def bam_file_template(self) -> Path:
        return Path(self.config["bam_file_template"]).resolve()

    @property
    def bai_file_template(self) -> Path:
        return self.bam_file_template.with_suffix(".bam.bai")

    @property
    def read_counts_file_template(self) -> Path:
        return Path(self.config["read_counts_file_template"]).resolve()

    @property
    def control_bam_file_template(self) -> Path:
        return Path(self.config["control_bam_file_template"]).resolve()

    @property
    def control_bai_file_template(self) -> Path:
        return self.control_bam_file_template.with_suffix(".bam.bai")

    @property
    def snp_file_template(self) -> Path:
        return Path(self.config["snp_file_template"]).resolve()

    @property
    def clone_filter_file_template(self) -> Path:
        return Path(self.config["clone_filter_file_template"]).resolve()

    @property
    def hapclone_results_file_template(self) -> Path:
        return Path(self.config["hapclone_results_file_template"]).resolve()

    # INPUT PREPROC DATA

    @property
    def black_list_file(self) -> Path:
        return Path(self.config["black_list_file"]).resolve()

    @property
    def centromere_file(self) -> Path:
        return Path(self.config["centromere_file"]).resolve()

    @property
    def map_file(self) -> Path:
        return Path(self.config["map_file"]).resolve()

    @property
    def ref_genome_file(self) -> Path:
        return Path(self.config["ref_genome_file"]).resolve()

    # DOWN SAMPLE BAM SETTINGS

    @property
    def read_length(self) -> int:
        """Used to compute number of reads needed to down sample bam file to"""
        return self.config["read_length"]

    @property
    def genome_length(self) -> int:
        """Used to compute number of reads needed to down sample bam file to"""
        return self.config["genome_length"]

    # PREPROC SETTINGS

    @property
    def bin_size(self):
        return self.config.get("bin_size", "500kb")

    @property
    def add_chr_prefix(self):
        return self.config.get("add_chr_prefix", True)

    @property
    def min_bqual(self):
        return int(self.config.get("min_bqual", 30))

    @property
    def min_mqual(self):
        return int(self.config.get("min_mqual", 30))

    @property
    def chromosomes(self):

        chroms = self.config.get("chromosomes", [x for x in range(1, 23)] + ["X"])

        chroms = [str(x) for x in chroms]

        if "autosomes" in self.config["chromosomes"]:
            
            chroms.remove("autosomes")
            
            chroms = [str(x) for x in range(1, 23)] + chroms

        for c in chroms:
            
            assert c in [str(x) for x in range(1, 23)] + ["X"]

        if self.add_chr_prefix:
            
            chroms = ["chr{}".format(x) for x in chroms]

        return chroms

    # CFCLONE SETTINGS

    @property
    def cfclone_config(self):
        return {
            "ctdna_file": str(self.cfclone_ctdna_template),
            "clone_cn_file": str(self.cfclone_clone_cn_template),
            "clone_tree_newick": str(None),  # cfclone-smk needs o.w. schema validation
            "num_chains": self.num_chains,
            "num_rounds": self.num_rounds,
            "num_threads": self.num_threads,
            "num_restarts": self.num_model_replicates,
            "out_dir": "<out_dir>",
            "pipeline_dir": "<pipeline_dir>",
        }

    @property
    def num_chains(self):
        return self.config.get("num_chains", 16)

    @property
    def num_rounds(self):
        return self.config.get("num_rounds", 10)

    @property
    def num_threads(self):
        return self.config.get("num_threads", 1)

    @property
    def cfclone_use_outlier(self):
        return self.config["cfclone_use_outlier"]

    # OUTPUTS

    @property
    def out_dir(self) -> Path:
        return Path(self.config["out_dir"]).resolve()

    @property
    def pipeline_dir(self) -> Path:
        return Path(self.config["pipeline_dir"]).resolve()

    @property
    def copied_config(self) -> Path:
        return self.out_dir.joinpath("config.yaml")

    # DOWN SAMPLING AND MIXING BAMS OUTPUTS

    @property
    def down_sampled_bam_file_template(self) -> Path:
        return self.pipeline_dir.joinpath(
            "coverage_{coverage_id}",
            # "final_sample_{final_sample_id}",
            "mixture_{mixture_id}",
            "proportion_{proportion_id}",
            "data_seed_{data_seed_id}",
            "{bam_id}.bam",
        )

    @property
    def mixed_bam_file_template(self) -> Path:
        return self.out_dir.joinpath(
            "coverage_{coverage_id}",
            # "final_sample_{final_sample_id}",
            "mixture_{mixture_id}",
            "proportion_{proportion_id}",
            "data_seed_{data_seed_id}",
            "mixed.bam",
        )

    @property
    def down_sampled_bai_file_template(self) -> Path:
        return self.down_sampled_bam_file_template.with_suffix(".bam.bai")

    @property
    def down_sampled_total_reads_template(self) -> Path:
        return self.down_sampled_bam_file_template.with_suffix(".tsv")

    @property
    def mixed_bai_file_template(self) -> Path:
        return self.mixed_bam_file_template.with_suffix(".bam.bai")

    @property
    def mixed_bam_total_reads_template(self) -> Path:
        return self.mixed_bam_file_template.with_suffix(".tsv")

    # PREPROC OUTPUTS

    @property
    def ref_dir(self) -> Path:
        return self.pipeline_dir.joinpath("ref")

    @property
    def gc_template(self):
        return self.ref_dir.joinpath("gc.tsv")

    @property
    def map_template(self):
        return self.ref_dir.joinpath("map.tsv")

    @property
    def working_dir(self) -> Path:
        return self.pipeline_dir.joinpath(
            "working",
            "coverage_{coverage_id}",
            # "final_sample_{final_sample_id}",
            "mixture_{mixture_id}",
            "proportion_{proportion_id}",
            "data_seed_{data_seed_id}",
        )

    @property
    def allele_counts_chrom_template(self):
        return self.working_dir.joinpath(
            "allele_counts_chrom", "{chrom}", "allele_counts.tsv.gz"
        )

    @property
    def reads_chrom_template(self):
        return self.working_dir.joinpath("reads_chrom", "{chrom}", "read_counts.tsv.gz")

    @property
    def normal_reads_chrom_template(self):
        return self.working_dir.joinpath(
            "reads_chrom", "normal", "{chrom}", "read_counts.tsv.gz"
        )

    @property
    def allele_counts_template(self):
        return self.working_dir.joinpath("allele_counts", "allele_counts.tsv.gz")

    @property
    def baf_template(self):
        return self.working_dir.joinpath("baf", "baf.tsv.gz")

    @property
    def rdr_template(self):
        return self.working_dir.joinpath("rdr", "rdr.tsv.gz")

    @property
    def reads_template(self):
        return self.working_dir.joinpath("reads", "read_counts.tsv.gz")

    @property
    def normal_reads_template(self):
        return self.working_dir.joinpath("reads", "normal", "read_counts.tsv.gz")

    @property
    def reads_wig_template(self):
        return self.working_dir.joinpath("reads", "read_counts.wig")

    @property
    def combined_results_template(self):
        return self.working_dir.joinpath("counts.tsv.gz")

    @property
    def baf_plot_template(self):
        return self.working_dir.joinpath("plots", "baf", "baf.png")

    @property
    def rdr_plot_template(self):
        return self.working_dir.joinpath("plots", "rdr", "rdr.png")

    # CFCLONE INPUT DATA

    @property
    def cfclone_input_dir(self):
        return self.out_dir.joinpath("input")

    @property
    def cfclone_clone_cn_template(self):
        return self.cfclone_input_dir.joinpath("clone_cn", "clone_cn.tsv.gz")

    @property
    def cfclone_ctdna_template(self):
        return self.cfclone_input_dir.joinpath(
            "ctdna",
            "coverage_{coverage_id}",
            # "final_sample_{final_sample_id}",
            "mixture_{mixture_id}",
            "proportion_{proportion_id}",
            "data_seed_{data_seed_id}",
            "data.tsv.gz",
        )
    
    # CFCLONE-SMK: INPUTS AND OUTPUTS

    @property
    def cfclone_out_dir(self):
        return self.out_dir.joinpath(
            "coverage_{coverage_id}",
            # "final_sample_{final_sample_id}",
            "mixture_{mixture_id}",
            "proportion_{proportion_id}",
            "data_seed_{data_seed_id}",
            "out_dir",
        )

    @property
    def cfclone_pipeline_dir(self):
        return self.out_dir.joinpath(
            "coverage_{coverage_id}",
            # "final_sample_{final_sample_id}",
            "mixture_{mixture_id}",
            "proportion_{proportion_id}",
            "data_seed_{data_seed_id}",
            "pipeline_dir",
        )

    @property
    def cfclone_clone_cn_template(self):
        return self.out_dir.joinpath(
            "input", "clone_cn", "clone_cn.tsv.gz"
        )

    @property
    def cfclone_ctdna_template(self):
        return self.out_dir.joinpath(
            "input",
            "ctdna",
            "coverage_{coverage_id}",
            # "final_sample_{final_sample_id}",
            "mixture_{mixture_id}",
            "proportion_{proportion_id}",
            "data_seed_{data_seed_id}",
            "data.tsv.gz",
        )

    @property
    def experiment_configuration(self):
        return self.cfclone_out_dir.joinpath("config.yaml")

    @property
    def merged_tumour_content_file(self):
        return self.cfclone_out_dir.joinpath("tumour_content.tsv")

    @property
    def merged_evidence_file(self):
        return self.cfclone_out_dir.joinpath("evidence.tsv")

    @property
    def merged_summary_file(self):
        return self.cfclone_out_dir.joinpath("summary.tsv")

    # SUMMARY OUTPUTS

    @property
    def outputs(self) -> Path:
        return self.out_dir.joinpath("outputs")

    @property
    def summary_file(self):
        return self.outputs.joinpath("summary.tsv")

    @property
    def down_sampled_summary_file(self) -> Path:
        return self.outputs.joinpath("down_sampled_bams_summary.tsv")

    @property
    def mixed_bams_summary_file(self) -> Path:
        return self.outputs.joinpath("mixed_bams_summary.tsv")

    # @property
    # def plot_summary_file(self) -> Path:
    #     return self.outputs.joinpath("tfs.png")

    # @property
    # def plot_summary_file_w_cohort(self) -> Path:
    #     return self.outputs.joinpath("tfs_w_cohort.png")

    @property
    def pipeline_files(self) -> list[str]:

        files = []

        files.append(self.copied_config)

        files.append(self.mixed_bams_summary_file)

        files.append(self.summary_file)

        # files.append(self.plot_summary_file)

        # files.append(self.plot_summary_file_w_cohort)

        for c in self.coverage_ids:
            
            # for s in self.final_sample_ids:
            for s in self.mixture_ids:

                for i in self.proportion_ids:
                    
                    for r in self.data_seed_ids:
                        
                        files.append(
                            str(self.experiment_configuration).format(
                                coverage_id=c,
                                mixture_id=s,
                                proportion_id=i,
                                data_seed_id=r,
                            )
                        )

                        files.append(
                            str(self.merged_tumour_content_file).format(
                                coverage_id=c,
                                mixture_id=s,
                                proportion_id=i,
                                data_seed_id=r,
                            )
                        )

                        files.append(
                            str(self.merged_evidence_file).format(
                                coverage_id=c,
                                mixture_id=s,
                                proportion_id=i,
                                data_seed_id=r,
                            )
                        )

                        files.append(
                            str(self.merged_summary_file).format(
                                coverage_id=c,
                                mixture_id=s,
                                proportion_id=i,
                                data_seed_id=r,
                            )
                        )

                        # files.append(
                        #     str(self.rdr_plot_template).format(coverage_id=c, final_sample_id = s, proportion_id=i, data_seed_id=r)
                        # )

                        # files.append(
                        #     str(self.baf_plot_template).format(coverage_id=c, final_sample_id = s, proportion_id=i, data_seed_id=r)
                        # )

        return files

    # HELPERS FOR RULES

    # def get_sample(self, final_sample_id: int, bam_id: str) -> str:

    #     if bam_id == "initial":

    #         return self.initial_sample

    #     else:

    #         return self.final_samples[final_sample_id]
        
    def get_bam_file(self, wildcards: dict) -> str:
        return str(self.bam_file_template).format(
            patient=self.patient, 
            sample=int(self.mixtures[wildcards.mixture_id][wildcards.bam_id]),
        )

    def get_bai_file(self, wildcards: dict) -> str:
        return str(Path(self.get_bam_file(wildcards)).with_suffix('.bam.bai'))

    def get_mixed_samples(self, wildcards: dict) -> str:
        return (
            self.initial_sample
            + "+"
            + self.final_samples[int(wildcards.final_sample_id)]
        )

    def get_bams_to_mix(self, wildcards: dict) -> list[str]:
        
        coverage_id = int(wildcards.coverage_id)
        
        mixture_id = int(wildcards.mixture_id)
        
        proportion_id = int(wildcards.proportion_id)
        
        # final_sample_id = int(wildcards.final_sample_id)
        
        data_seed_id = int(wildcards.data_seed_id)

        files = []

        for b in self.bam_ids:
            
            p = self.proportions[proportion_id] if b == 'final' else 1 - p

            if p > 0.:

                files.append(
                    str(self.down_sampled_bam_file_template).format(
                        coverage_id=coverage_id,
                        # final_sample_id=final_sample_id,
                        mixture_id=mixture_id,
                        proportion_id=proportion_id,
                        data_seed_id=data_seed_id,
                        bam_id=b,
                    )
                )

        return files

    @property
    def gather_down_sampled_summary_files(self) -> list[str]:

        files = []

        for c in self.coverage_ids:

            for p in self.proportion_ids:
                
                for r in self.data_seed_ids:

                    for b in self.bam_ids:

                        p = self.compute_bam_proportion(
                            bam_id=b,
                            coverage_id=c,
                            proportion_id=p,
                        )

                        if p > 0.0:

                            files.append(
                                str(self.down_sampled_total_reads_template).format(
                                    coverage_id=c, 
                                    proportion_id=p,
                                    data_seed_id=r,
                                    bam_id=b,
                                )
                            )

        return files

    @property
    def get_read_counts_file(self) -> str:
        return str(self.read_counts_file_template).format(patient=self.patient)

    @property
    def get_control_bam_file(self) -> str:
        return str(self.control_bam_file_template).format(patient=self.patient)

    @property
    def get_control_bai_file(self) -> str:
        return str(self.control_bai_file_template).format(patient=self.patient)

    @property
    def get_snp_file(self) -> str:
        return str(self.snp_file_template).format(patient=self.patient)

    @property
    def get_hapclone_results_file(self) -> str:
        return str(self.hapclone_results_file_template).format(patient=self.patient)

    @property
    def get_clone_filter_file(self) -> str:
        return str(self.clone_filter_file_template).format(patient=self.patient)

    # HELPERS FOR PREPROC RULES

    @staticmethod
    def parse_region_size(size_str):
        return int(size_str.replace("kb", "")) * int(1e3)

    # HELPERS TO GATHERS FILES

    @property
    def gather_mixed_bams_summary_files(self) -> list[str]:

        files = []

        for c in self.coverage_ids:

            # for s in self.final_sample_ids:
            
            for s in self.mixture_ids:

                for i in self.proportion_ids:
                    
                    for r in self.data_seed_ids:

                        file = str(self.mixed_bam_total_reads_template).format(
                            coverage_id=c,
                            # final_sample_id=s,
                            mixture_id=s,
                            proportion_id=i,
                            data_seed_id=r,
                        )

                        files.append(file)

        return files

    def gather_files_by_chrom(self, file_template: Path, wildcards: dict) -> list[str]:

        files = []

        for c in self.chromosomes:

            files.append(
                str(file_template).format(
                    coverage_id=int(wildcards["coverage_id"]),
                    mixture_id=int(wildcards['mixture_id']),
                    proportion_id=int(wildcards["proportion_id"]),
                    data_seed_id=int(wildcards["data_seed_id"]),
                    chrom=c,
                )
            )

        return files

    # HELPERS TO COMPUTE FRACTION OF BAM FILE NEED TO KEEP

    def compute_down_sample_proportion(self, wildcards: dict) -> str:

        bam_file_prop = self.compute_bam_proportion(
            bam_id=wildcards.bam_id,
            coverage_id=int(wildcards.coverage_id),
            mixture_id=int(wildcards.mixture_id),
            proportion_id=int(wildcards.proportion_id),
        )

        return format(bam_file_prop, "f")  # store as string for cli

    def compute_bam_proportion(
        self,
        bam_id: str,
        coverage_id: int,
        mixture_id: int,
        proportion_id: int
    ) -> float:

        # LOAD TOTAL NUMBER OF READS IN BAM FILE

        sample_name = self.mixtures[mixture_id][bam_id]
        
        num_reads_avail = (
            
            pd.read_csv(self.get_read_counts_file, sep="\t")
            
            .loc[lambda df: df["file"] == "{}.bam".format(sample_name), "reads"]
            
            .values[0]
            
        )

        # COMPUTE TOTAL NUMBER OF READS NEEDED FROM BAM FILE

        cov = self.coverages[coverage_id]
        
        prop = self.proportions[proportion_id]
        
        p = prop if bam_id == 'final' else 1 - prop

        num_reads_needed = int((self.genome_length / self.read_length) * cov * p)

        # CHECK IF THERE ARE ENOUGH READS IN THE BAM FILE

        if num_reads_avail < num_reads_needed:

            msg = "Reads needed is greater than reads avaible: {n} > {a}"

            raise ValueError(msg.format(n=num_reads_needed, a=num_reads_avail))

        # PROPORTION OF BAM FILE WE DOWN SAMPLE TO

        bam_file_prop = num_reads_needed / num_reads_avail

        return bam_file_prop

    # HELPERS FOR LOG AND BENCHMARK FILES

    @property
    def log_dir(self) -> Path:
        return self.out_dir.joinpath("logs")

    @property
    def benchmark_dir(self) -> Path:
        return self.out_dir.joinpath("benchmarks")

    def get_log_file(self, template: Path) -> Path:
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.log_dir.joinpath(parent, rel_path)

    def get_benchmark_file(self, template: Path) -> Path:
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.benchmark_dir.joinpath(parent, rel_path)

    def _get_relative_path(self, template: Path) -> tuple[str, Path]:
        try:
            rel_path = template.relative_to(self.pipeline_dir)
            parent = "pipeline"
        except ValueError:
            rel_path = template.relative_to(self.out_dir)
            parent = "output"
        return parent, rel_path

    # HELPERS FOR EMAIL NOTIFICATIONS

    @property
    def email(self) -> str:
        return self.config.get("email", "lepurmatteo@gmail.com")

    def email_notification(
        self, on: str, workflow: str, configfile: str, imgs: list[str] | None = None
    ) -> None:

        configfile = Path(configfile).resolve()

        msg_template = "configfile:{config}"

        msg = msg_template.format(config=configfile)

        subj_template = "-s {wf}:{n}"

        subj = subj_template.format(wf=workflow, n=on)

        if (on == "success") and (imgs is not None):

            img_template = "-a {img} "

            imgs = self.get_imgs_to_send(imgs)

            att = ""

            for i in imgs:

                att0 = img_template.format(img=i)

                att += att0

            cmd_template = "echo {msg} | mail {sub} {att} {email}"

            cmd = cmd_template.format(msg=msg, sub=subj, att=att, email=self.email)

        else:

            cmd_template = "echo {msg} | mail {sub} {email}"

            cmd = cmd_template.format(msg=msg, sub=subj, email=self.email)

        shell(cmd)

    @staticmethod
    def get_imgs_to_send(imgs: list[str]) -> list[str]:

        LIMIT = 10240000

        ALLOWED = LIMIT * 0.66

        total_size = 0

        allowed_imgs = []

        for img in imgs:

            file_path = Path(img)

            file_size = file_path.stat().st_size

            total_size += file_size

            if total_size < ALLOWED:

                allowed_imgs.append(file_path)

            else:

                total_size -= file_size

        return allowed_imgs
