from janis_bioinformatics.tools import BioinformaticsWorkflow
from janis_core import String, Array, File, ScatterDescription, ScatterMethod

from janis_bioinformatics.data_types import FastqGzPair, FastaWithDict

from janis_pipelines.alignment import AgilentUMIalignment


class AgilentUMIalignmentMultiSample(BioinformaticsWorkflow):
    def id(self):
        return "AgilentUMIalignmentMultiSample"

    def friendly_name(self):
        return "Preprocessing of Agilent UMI data from FASTQ to UMI duplicate marked bam for multiple samples"

    def version(self):
        return "0.0.0"

    def constructor(self):
        self.input("reads", Array(Array(FastqGzPair)))
        self.input("sample_name", Array(String()))

        self.input("reference", FastaWithDict())

        self.input("referenceAlt", File())

        # For Agent Trimmer
        self.input("agentlibrary", String())
        # For CutAdapt
        self.input("cutadapt_adapters", File(optional=True))

        ## STEPS
        self.step(
            "alignment_workflow",
            AgilentUMIalignment(
                reads=self.reads,
                sample_name=self.sample_name,
                reference=self.reference,
                referenceAlt=self.referenceAlt,
                agentlibrary=self.agentlibrary,
                cutadapt_adapters=self.cutadapt_adapters,
            ),
            scatter=ScatterDescription(
                ["reads", "sample_name"],
                method=ScatterMethod.dot,
                labels=self.sample_name,
            ),
        )

        ## OUTPUTs
        self.capture_outputs_from_step(self.alignment_workflow)

