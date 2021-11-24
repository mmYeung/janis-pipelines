from typing import Optional, List
from janis_core import String, Array, File, InputDocumentation

from janis_bioinformatics.data_types import FastqGzPair, FastaWithDict

from janis_bioinformatics.tools.babrahambioinformatics import FastQC_0_11_8

from janis_bioinformatics.tools.pmac import ParseFastqcAdaptors

from janis_bioinformatics.tools.bwakit import BwaPostAltAligner

from janis_bioinformatics.tools.common import MergeAndMarkBamsUMI_4_1_2


from janis_pipelines import CaptureSomaticTumourOnlyMultiCallersVariantsOnly


class CaptureSomaticTumourOnly(
    CaptureSomaticTumourOnlyMultiCallersVariantsOnly
):
    def id(self):
        return "CaptureSomaticTumourOnly"

    def friendly_name(self):
        return "Capture Somatic Tumour Only (Multi callers)"

    def version(self):
        return "0.0.0"

    def constructor(self):
        self.add_inputs()

        self.step("fastqc", FastQC_0_11_8(reads=self.reads, scatter="reads"))

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
        super.add_inputs()
        self.input("reads", Array(FastqGzPair))
        self.input("referenceAlt", File())

    def add_preprocessing(self):

        sub_inputs = {
            "reference": self.reference,
            "referenceAlt": self.referenceAlt,
            "cutadapt_adapter": self.getfastqc_adapters,
            "cutadapt_removeMiddleAdapter": self.getfastqc_adapters,
        }

        self.step(
            "alignUMI",
            BwaPostAltAligner(
                fastqs=self.reads, sample_name=self.sample_name, **sub_inputs
            ),
            scatter=[
                "fastqs",
                "cutadapt_adapter",
                "cutadapt_removeMiddleAdapter",
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
