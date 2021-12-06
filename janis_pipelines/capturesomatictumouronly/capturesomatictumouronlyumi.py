from typing import Optional, List
from janis_core import String, Array, File, InputDocumentation

from janis_bioinformatics.data_types import FastqGzPair, FastaWithDict

from janis_bioinformatics.tools.agent import AgentTrimmer_2_0_2

from janis_bioinformatics.tools.babrahambioinformatics import FastQC_0_11_8

from janis_bioinformatics.tools.pmac import ParseFastqcAdaptors

from janis_bioinformatics.tools.bwakit import BwaPostAltAlignerUMI

from janis_bioinformatics.tools.common import MergeAndMarkBamsUMI_4_1_2


from janis_pipelines.capturesomatictumouronly.capturesomatictumouronly_variantsonly import (
    CaptureSomaticTumourOnlyMultiCallersVariantsOnly,
)


class CaptureSomaticTumourOnlyUMI(CaptureSomaticTumourOnlyMultiCallersVariantsOnly):
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
            AgentTrimmer_2_0_2(
                read1=self.reads[0],
                read2=self.reads[1],
                outdir="/tmp",
                library=self.agentlibrary,
                agentVersion="2.0.2",
            ),
            scatter=["read1", "read2"],
        )

        self.step("fastqc", FastQC_0_11_8(reads=self.agenttrim.out), scatter="reads")

        self.step(
            "getfastqc_adapters",
            ParseFastqcAdaptors(
                fastqc_datafiles=self.fastqc.datafile,
                cutadapt_adaptors_lookup=self.cutadapt_adapters,
            ),
            scatter="fastqc_datafiles",
        )

        self.add_preprocessing()

        # self.step(
        #     "callers",
        #     CaptureSomaticTumourOnlyMultiCallersVariantsOnly(
        #         bam=self.add_preprocessing.out_bam,
        #         referenceFolder=self.referenceFolder,
        #     ),
        # )

    def add_inputs(self):
        self.input("reads", Array(FastqGzPair))
        self.input("sample_name", String())
        self.input("referenceAlt", File())
        # For Agent Trimme
        self.input("agentlibrary", String())
        # For CutAdapt
        self.input("cutadapt_adapters", File(optional=True))
        self.add_inputs_for_reference()
        self.add_inputs_for_configuration()
        self.add_inputs_for_intervals()
        self.add_inputs_for_vc()

    def add_preprocessing(self):

        sub_inputs = {
            "reference": self.reference,
            "referenceAlt": self.referenceAlt,
            "cutadapt_adapter": self.getfastqc_adapters,
            "cutadapt_removeMiddle3Adapter": self.getfastqc_adapters,
        }

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

        self.output("out_bam", source=self.merge_and_mark.out)
        self.output(
            "umimetrics",
            source=self.merge_and_mark.umimetrics,
            output_folder=["stats"],
            output_name=self.sample_name + "umimetrics.txt",
        )
        self.output(
            "metrics",
            source=self.merge_and_mark.metrics,
            output_folder=["stats"],
            output_name=self.sample_name + "metrics.txt",
        )
        self.output("out_fastqc_reports", source=self.fastqc.out)
