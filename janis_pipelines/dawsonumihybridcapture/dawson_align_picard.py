from janis_bioinformatics.tools import BioinformaticsWorkflow

from janis_bioinformatics.data_types import (
    FastqGzPair,
    Bam,
    Vcf,
    BamBai,
    FastaWithDict,
)
from janis_bioinformatics.tools.babrahambioinformatics import FastQC_0_11_8
from janis_bionformatics.tools.pmac import ParseFastqcAdaptors
from janis_bioinformatics.tools.bwakit import BwaPostAltAlignerUMI
from janis_bioinformatics.tools.common.mergeandmarkumi import (
    MergeAndMarkBamsUMI_4_1_2,
)

from janis_core import String, Array, File


class AgilentSureSelectXTHS2(BioinformaticsWorkflow):
    def id(self):
        return "AgilentSureSelectXTHS2"

    def friendly_name(self):
        return "Agilent SureSelect XTHS2"

    def constructor(self):
        ## Inputs
        self.input("fastqs", Array(FastqGzPair), doc="Fastq Files")
        self.input("reference", FastaWithDict)
        self.input("")

        ## Steps
        self.step("fastqc", FastQC_0_11_8(reads=self.fastqs), scatter="reads")

        self.step(
            "getfastqc_adapters",
            ParseFastqcAdaptors(
                fastqc_datafiles=self.fastqc.datafile,
                cutadapt_adaptors_lookup=self.cutadapt_adapters,
            ),
            scatter="fastqc_datafiles",
        )

        self.step(
            "align",
            BwaPostAltAlignerUMI(
                fastq=self.fastqs,
                reference=self.reference,
                referenceAlt=self.referenceAlt,
                sample_name=self.sample_name,
                sortsam_tmpDir="./tmp",
                cutadapt_adapter=self.getfastqc_adapters,
                cutadapt_rmoveMiddle3Adapter=self.getfastqc_adapters,
            ),
            scatter=[
                "fastq",
                "cutadapt_adapter",
                "cutadapt_removeMiddle3Adapter",
            ],
        )

        self.step(
            "merge_fixmateinfo_mark",
            MergeAndMarkBamsUMI_4_1_2(
                bams=self.align.out, sampleName=self.sample_name
            ),
        )

        self.output(
            "out_bam",
            source=self.merge_and_mark.out,
            output_folder="Bam",
            doc="Aligned and index bam.",
            output_name=self.sample_name,
        )

        ## Outputs
        self.step(
            "out_fastqc_reports",
            source=self.fastqc.out,
            output_folder="reports",
            doc="A zip file of the FastQC quality report.",
        )
