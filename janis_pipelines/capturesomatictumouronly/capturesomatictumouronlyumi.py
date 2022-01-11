from typing import Optional, List
from janis_bioinformatics.data_types.fastq import FastqPair
from janis_core import String, Array, File, WorkflowBuilder, StringFormatter

from janis_bioinformatics.data_types import FastqGzPair, FastaWithDict

from janis_bioinformatics.tools.agent import AgentTrimmer_2_0_2

from janis_bioinformatics.tools.babrahambioinformatics import FastQC_0_11_8

from janis_bioinformatics.tools.pmac import ParseFastqcAdaptors

from janis_bioinformatics.tools.bwakit import BwaPostAltAlignerUMI

from janis_bioinformatics.tools.common import MergeAndMarkBamsUMI_4_1_2

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
            "agenttrim",
            self.umi_trimmer_subworkflow(
                fastqPair=self.reads, agentlibrary=self.agentlibrary
            ),
            scatter=["fastqPair"],
        )

        self.step(
            "fastqc", FastQC_0_11_8(reads=self.agenttrim.out), scatter="reads"
        )

        self.step(
            "getfastqc_adapters",
            ParseFastqcAdaptors(
                fastqc_datafiles=self.fastqc.datafile,
                cutadapt_adaptors_lookup=self.cutadapt_adapters,
            ),
            scatter="fastqc_datafiles",
        )

        self.add_preprocessing()

        self.step(
            "callers",
            CaptureSomaticTumourOnlyMultiCallersVariantsOnly(
                bam=self.merge_and_mark.out,
                sample_name=self.sample_name,
                reference=self.reference,
                referenceFolder=self.referenceFolder,
                intervals=self.intervals,
                minVaf=self.minVaf,
                minMQ=self.minMQ,
                minAD=self.minAD,
                minDP=self.minDP,
                minBQ=self.minBQ,
                pileupMaxDepth=self.pileupMaxDepth,
                pileupMinBQ=self.pileupMinBQ,
                varscanPval=self.varscanPval,
                piscesVCminVQ=self.piscesVCminVQ,
                piscesVQRminVQ=self.piscesVQRminVQ,
                snps_dbsnp=self.snps_dbsnp,
                snps_1000gp=self.snps_1000gp,
                known_indels=self.known_indels,
                mills_indels=self.mills_indels,
                gnomad=self.gnomad,
                panel_of_normals=self.panel_of_normals,
            ),
        )

        ## Outputs
        self.output(
            "combined",
            source=self.callers.out_variants,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}.combined", self.sample_name
            ),
        )

        self.output(
            "vardict",
            source=self.callers.out_variants_vardict,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}.vardict.recode", samplename=self.sample_name
            ),
        )
        self.output(
            "varscan",
            source=self.callers.out_variants_varscan2,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}.varscan2.recode", samplename=self.sample_name
            ),
        )
        self.output(
            "mutect2",
            source=self.callers.out_variants_pass_gatk,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}.mutect2.recode", samplename=self.sample_name
            ),
        )
        self.output(
            "pisces",
            source=self.callers.out_variants_pisces,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}.pisces.recode", samplename=self.sample_name
            ),
        )

        self.output(
            "vardict_raw",
            source=self.callers.variants_vardict,
            output_folder=["variants", "raw"],
            output_name=StringFormatter(
                "{samplename}.vardict", samplename=self.sample_name
            ),
        )
        self.output(
            "varscan2_raw",
            source=self.callers.variants_varscan2,
            output_folder=["variants", "raw"],
            output_name=StringFormatter(
                "{samplename}.varscan", samplename=self.sample_name
            ),
        )
        self.output(
            "mutect2_raw",
            source=self.callers.out_variants_gatk,
            output_folder=["variants", "raw"],
            output_name=StringFormatter(
                "{samplename}.mutect2", samplename=self.sample_name
            ),
        )
        self.output(
            "pisces_raw",
            source=self.callers.variants_pisces,
            output_folder=["variants", "raw"],
            output_name=StringFormatter(
                "{samplename}.pisces", samplename=self.sample_name
            ),
        )

        self.output(
            "pisces_bam",
            source=self.callers.pisces_bam,
            output_folder=["Bam"],
            output_name=StringFormatter(
                "{samplename}.hygea.stitcher", samplename=self.sample_name
            ),
        )
        self.output(
            "gatk_bam",
            source=self.callers.gatk_bam,
            output_folder=["Bam"],
            output_name=StringFormatter(
                "{samplename}.gatk", samplename=self.sample_name
            ),
        )

    def add_inputs(self):
        self.input("reads", Array(FastqGzPair))
        self.input("sample_name", String())
        self.input("referenceAlt", File())
        # For Agent Trimmer
        self.input("agentlibrary", String())
        # For CutAdapt
        self.input("cutadapt_adapters", File(optional=True))
        self.add_inputs_for_reference()
        self.add_inputs_for_configuration()
        self.add_inputs_for_intervals()
        self.add_inputs_for_vc()

    @staticmethod
    def umi_trimmer_subworkflow(**connections):
        w = WorkflowBuilder("umi_trimmer_subworkflow")

        ## Inputs
        w.input("fastqPair", FastqGzPair())
        w.input("agentlibrary", String())

        w.step(
            "agenttrimsub",
            AgentTrimmer_2_0_2(
                read1=w.fastqPair[0],
                read2=w.fastqPair[1],
                outdir=".",
                library=w.agentlibrary,
                agentVersion="2.0.2",
            ),
        )

        w.output("out", source=w.agenttrimsub.out)

        return w(**connections)

    def add_preprocessing(self):

        sub_inputs = {
            "reference": self.reference,
            "referenceAlt": self.referenceAlt,
            "cutadapt_adapter": self.getfastqc_adapters,
            "cutadapt_removeMiddle3Adapter": self.getfastqc_adapters,
        }

        ## STEPS
        self.step(
            "alignUMI",
            BwaPostAltAlignerUMI(
                fastq=self.agenttrim.out,
                sample_name=self.sample_name,
                addumis=True,
                **sub_inputs,
            ),
            scatter=[
                "fastq",
                "cutadapt_adapter",
                "cutadapt_removeMiddle3Adapter",
            ],
        )

        self.step(
            "merge_and_mark",
            MergeAndMarkBamsUMI_4_1_2(
                bams=self.alignUMI.out, sample_name=self.sample_name
            ),
        )

        ## OUTPUTS
        self.output(
            "out_bam",
            source=self.merge_and_mark.out,
            output_folder=["Bam"],
            output_name=StringFormatter("{samplename}.umarkdups"),
        )
        self.output(
            "umimetrics",
            source=self.merge_and_mark.umimetrics,
            output_folder=["stats"],
            output_name=self.sample_name + "umimetrics",
        )
        self.output(
            "metrics",
            source=self.merge_and_mark.metrics,
            output_folder=["stats"],
            output_name=self.sample_name + "metrics",
        )
        self.output("out_fastqc_reports", source=self.fastqc.out)
