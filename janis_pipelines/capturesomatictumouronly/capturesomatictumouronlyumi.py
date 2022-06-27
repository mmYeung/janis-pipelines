from typing import Optional, List
from janis_bioinformatics.data_types.fastq import FastqPair
from janis_core import String, Array, File, WorkflowBuilder, StringFormatter

from janis_bioinformatics.data_types import FastqGzPair, FastaWithDict

from janis_pipelines.alignment.agilentUMIalignment import AgilentUMIalignment

from janis_pipelines.capturesomatictumouronly.capturesomatictumouronly_variantsonly import (
    CaptureSomaticTumourOnlyMultiCallersVariantsOnly,
)


class CaptureSomaticTumourOnlyUMI(
    CaptureSomaticTumourOnlyMultiCallersVariantsOnly
):
    def id(self):
        return "CaptureSomaticTumourOnlyUMI"

    def friendly_name(self):
        return "Capture Somatic Tumour Only with UMIs (Multi callers)"

    def version(self):
        return "0.0.0"

    def constructor(self):
        self.add_inputs()

        self.step(
            "preprocessing",
            AgilentUMIalignment(
                reads=self.reads,
                sample_name=self.sample_name,
                reference=self.reference,
                reference_alt=self.reference_alt,
                agentlibrary=self.agentlibrary,
                cutadapt_adapters=self.cutadapt_adapters,
            ),
        )

        self.step(
            "callers",
            CaptureSomaticTumourOnlyMultiCallersVariantsOnly(
                bam=self.preprocessing.out_bam,
                sample_name=self.sample_name,
                reference=self.reference,
                reference_folder=self.reference_folder,
                snps_dbsnp=self.snps_dbsnp,
                snps_1000gp=self.snps_1000gp,
                known_indels=self.known_indels,
                mills_indels=self.mills_indels,
                gnomad=self.gnomad,
                panel_of_normals=self.panel_of_normals,
                intervals=self.intervals,
                min_vaf=self.min_vaf,
                min_mq=self.min_mq,
                min_ad=self.min_ad,
                min_dp=self.min_dp,
                min_bq=self.min_bq,
                pileup_max_depth=self.pileup_max_depth,
                pileup_min_bq=self.pileup_min_bq,
                varscan_pval=self.varscan_pval,
                pisces_vc_min_vq=self.pisces_vc_min_vq,
                pisces_vqr_min_vq=self.pisces_vqr_min_vq,
                pisces_awk_script=self.pisces_awk_script,
            ),
        )

        ## Outputs
        self.output(
            "umarkdups_cram",
            source=self.preprocessing.out_cram,
            output_folder=["Cram"],
            output_name=StringFormatter(
                "{samplename}_umarkdups", samplename=self.sample_name
            ),
        )
        self.output(
            "umimetrics",
            source=self.preprocessing.umimetrics,
            output_folder=["stats"],
            output_name=StringFormatter(
                "{samplename}_umimetrics", samplename=self.sample_name
            ),
        )
        self.output(
            "metrics",
            source=self.preprocessing.metrics,
            output_folder=["stats"],
            output_name=StringFormatter(
                "{samplename}_metrics", samplename=self.sample_name
            ),
        )
        self.output(
            "out_fastqc_reports",
            source=self.preprocessing.out_fastqc_reports,
            output_folder=["QC"],
            output_name=StringFormatter(
                "{samplename}_fastqc_report", samplename=self.sample_name
            ),
        )

        self.output(
            "combined",
            source=self.callers.out_variants,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}_combined.vcf", samplename=self.sample_name
            ),
        )

        self.output(
            "vardict",
            source=self.callers.out_variants_vardict,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}_vardict.recode", samplename=self.sample_name
            ),
        )
        self.output(
            "varscan",
            source=self.callers.out_variants_varscan2,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}_varscan2.recode", samplename=self.sample_name
            ),
        )
        self.output(
            "mutect2",
            source=self.callers.out_variants_gatk,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}_mutect2.recode", samplename=self.sample_name
            ),
        )
        self.output(
            "pisces",
            source=self.callers.out_variants_pisces,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}_pisces.recode", samplename=self.sample_name
            ),
        )

        self.output(
            "vardict_raw",
            source=self.callers.variants_vardict,
            output_folder=["variants", "raw"],
            output_name=StringFormatter(
                "{samplename}_vardict", samplename=self.sample_name
            ),
        )
        self.output(
            "varscan2_raw",
            source=self.callers.variants_varscan2,
            output_folder=["variants", "raw"],
            output_name=StringFormatter(
                "{samplename}_varscan", samplename=self.sample_name
            ),
        )
        self.output(
            "mutect2_raw",
            source=self.callers.variants_gatk,
            output_folder=["variants", "raw"],
            output_name=StringFormatter(
                "{samplename}_mutect2", samplename=self.sample_name
            ),
        )
        self.output(
            "pisces_raw",
            source=self.callers.variants_pisces,
            output_folder=["variants", "raw"],
            output_name=StringFormatter(
                "{samplename}_pisces", samplename=self.sample_name
            ),
        )

        self.output(
            "pisces_cram",
            source=self.callers.pisces_cram,
            output_folder=["Cram"],
            output_name=StringFormatter(
                "{samplename}_hygea.stitcher", samplename=self.sample_name
            ),
        )
        self.output(
            "gatk_cram",
            source=self.callers.gatk_cram,
            output_folder=["Cram"],
            output_name=StringFormatter(
                "{samplename}_gatk", samplename=self.sample_name
            ),
        )

    def add_inputs(self):
        self.input("reads", Array(FastqGzPair))
        self.input("sample_name", String())
        self.input("reference_alt", File())
        # For Agent Trimmer
        self.input("agentlibrary", String())
        # For CutAdapt
        self.input("cutadapt_adapters", File(optional=True))
        self.add_inputs_for_reference()
        self.add_inputs_for_configuration()
        self.add_inputs_for_intervals()
        self.add_inputs_for_vc()

