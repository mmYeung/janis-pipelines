from janis_bioinformatics.tools import BioinformaticsWorkflow
from janis_core import String, Array, File, StringFormatter, WorkflowBuilder

from janis_bioinformatics.data_types import FastqGzPair, FastaWithDict

from janis_bioinformatics.tools.agent import AgentTrimmer_2_0_2

from janis_bioinformatics.tools.babrahambioinformatics import FastQC_0_11_8

from janis_bioinformatics.tools.pmac import ParseFastqcAdaptors

from janis_bioinformatics.tools.bwakit import BwaPostAltAlignerUMI

from janis_bioinformatics.tools.common import MergeAndMarkBamsUMI_4_1_2


class AgilentUMIalignment(BioinformaticsWorkflow):
    def id(self):
        return "AgilentUMIalignment"

    def friendly_name(self):
        return "Preprocessing of Agilent UMI data from FASTQ to UMI duplicate marked bam"

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

    def add_inputs(self):
        self.input("reads", Array(FastqGzPair))
        self.input("sample_name", String())

        self.input("reference", FastaWithDict)

        self.input("referenceAlt", File())

        # For Agent Trimmer
        self.input("agentlibrary", String())
        # For CutAdapt
        self.input("cutadapt_adapters", File(optional=True))

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
            output_name=StringFormatter(
                "{samplename}.umarkdups", samplename=self.sample_name
            ),
        )
        self.output(
            "umimetrics",
            source=self.merge_and_mark.umimetrics,
            output_folder=["stats"],
            output_name=StringFormatter(
                "{samplename}.umimetrics", samplename=self.sample_name
            ),
        )
        self.output(
            "metrics",
            source=self.merge_and_mark.metrics,
            output_folder=["stats"],
            output_name=StringFormatter(
                "{samplename}.metrics", samplename=self.sample_name
            ),
        )
        self.output(
            "out_fastqc_reports",
            source=self.fastqc.out,
            output_folder=["QC"],
            output_name=StringFormatter(
                "{samplename}_fastqc_report", samplename=self.sample_name
            ),
        )
