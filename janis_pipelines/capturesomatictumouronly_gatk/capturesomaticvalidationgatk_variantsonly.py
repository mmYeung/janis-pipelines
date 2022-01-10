from janis_core import (
    String,
    Array,
    Boolean,
    WorkflowMetadata,
    StringFormatter,
    InputQualityType,
)
from janis_unix.tools import UncompressArchive

from janis_bioinformatics.data_types import (
    BamBai,
    Bed,
    Vcf,
    CompressedVcf,
    VcfTabix,
)

from janis_pipelines.wgs_somatic_gatk.wgssomaticgatk_variantsonly import (
    WGSSomaticGATKVariantsOnly,
)

from janis_pipelines.reference import WGS_INPUTS

INPUT_DOCS = {
    **WGS_INPUTS,
    "normal_inputs": {
        "doc": "An array of NORMAL FastqGz pairs. These are aligned separately and merged "
        "to create higher depth coverages from multiple sets of reads",
        "quality": InputQualityType.user,
        "example": [
            ["normal_R1.fastq.gz", "normal_R2.fastq.gz"],
            ["normal_R1-TOPUP.fastq.gz", "normal_R2-TOPUP.fastq.gz"],
        ],
    },
    "tumor_inputs": {
        "doc": "An array of TUMOR FastqGz pairs. These are aligned separately and merged "
        "to create higher depth coverages from multiple sets of reads",
        "quality": InputQualityType.user,
        "example": [
            ["tumor_R1.fastq.gz", "tumor_R2.fastq.gz"],
            ["tumor_R1-TOPUP.fastq.gz", "tumor_R2-TOPUP.fastq.gz"],
        ],
    },
    "normal_name": {
        "doc": "Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem",
        "quality": InputQualityType.user,
        "example": "NA12878_normal",
    },
    "tumor_name": {
        "doc": "Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem",
        "quality": InputQualityType.user,
        "example": "NA12878_tumor",
    },
    "normal_bam": {
        "doc": "Indexed NORMAL bam to call somatic variants against",
        "quality": InputQualityType.user,
        "example": "NA12878-normal.bam",
    },
    "tumor_bam": {
        "doc": "Indexed TUMOR bam to call somatic variants against",
        "quality": InputQualityType.user,
        "example": "NA12878-normal.bam",
    },
}

from janis_bioinformatics.tools.dawson import GenerateChromosomeIntervalsFromBed

from janis_bioinformatics.tools.common.gatkbasecalbam import (
    GATKBaseRecalBQSRWorkflow_4_1_3,
)
from janis_bioinformatics.tools.variantcallers import (
    GatkSomaticVariantCallerTumorOnlyTargeted,
)
from janis_bioinformatics.tools.bcftools import (
    BcfToolsConcat_1_9,
    BcfToolsSort_1_9,
    BcfToolsNorm_1_9,
)

from janis_bioinformatics.tools.vcftools import VcfToolsvcftools_0_1_16


class CaptureSomaticValidationGATKVariantsOnly(WGSSomaticGATKVariantsOnly):
    def id(self):
        return "CaptureSomaticValidationGATKVariantsOnly"

    def friendly_name(self):
        return "Capture Somatic Validation(GATK only) [VARIANTS only]"

    def version(self):
        return "1.4.0"

    def add_inputs_for_intervals(self):
        self.input("intervals", Bed())

    def add_inputs(self):
        self.input("bam", BamBai())
        self.input("sample_name", String())

        self.add_inputs_for_configuration()
        self.add_inputs_for_intervals()
        self.add_inputs_for_reference()

    def add_inputs_for_configuration(self):
        self.input("gnomad", VcfTabix(), doc=INPUT_DOCS["gnomad"])
        self.input(
            "panel_of_normals", VcfTabix(), doc=INPUT_DOCS["panel_of_normals"]
        )
        self.input("genotype_germline", Boolean())

    def add_gatk_variantcaller(self, tumour_bam_source):

        recal_ins = {
            "reference": self.reference,
            "intervals": self.intervals,
            "snps_dbsnp": self.snps_dbsnp,
            "snps_1000gp": self.snps_1000gp,
            "known_indels": self.known_indels,
            "mills_indels": self.mills_indels,
        }

        self.step(
            "bqsr_tumour",
            GATKBaseRecalBQSRWorkflow_4_1_3(bam=tumour_bam_source, **recal_ins),
        )

        self.step(
            "vc_gatk",
            GatkSomaticVariantCallerTumorOnlyTargeted(
                bam=self.bqsr_tumour.out,
                intervals=self.intervals,
                reference=self.reference,
                gnomad=self.gnomad,
                genotype_germline=self.genotype_germline,
            ),
        )

        self.step("vc_gatk_sort", BcfToolsSort_1_9(vcf=self.vc_gatk.variants))

        self.step(
            "vc_gatk_normalise", BcfToolsNorm_1_9(vcf=self.vc_gatk_sort.out)
        )

        self.step(
            "vc_gatk_uncompressvcf",
            UncompressArchive(file=self.vc_gatk_normalise.out),
        )

        self.step(
            "vc_gatk_filterpass",
            VcfToolsvcftools_0_1_16(
                vcf=self.vc_gatk_uncompressvcf.out.as_type(Vcf),
                removeFileteredAll=True,
                recode=True,
                recodeINFOAll=True,
            ),
        )

    def constructor(self):
        ## INPUTS
        self.add_inputs()

        ## STEPS
        self.add_gatk_variantcaller(tumour_bam_source=self.bam)

        ## OUTPUTS
        self.output(
            "out_variants_gatk",
            source=self.vc_gatk_sort.out,
            output_folder=["variants", "unfiltered"],
            doc="Merged variants from the GATK caller",
        )

        self.output(
            "out_variants_pass_gatk",
            source=self.vc_gatk_filterpass.out,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}.gatk.recode.vcf", samplename=self.sample_name
            ),
        )

        self.output(
            "gatk_bam",
            source=self.vc_gatk.out_bam,
            output_folder=["Bam"],
            output_name=StringFormatter(
                "{samplename}.gatk.bam", samplename=self.sample_name
            ),
        )

    def bind_metadata(self):
        from datetime import date

        return WorkflowMetadata(
            version="1.0.0",
            contributors=["Miriam M Yeung"],
            dateCreated=date(2021, 11, 5),
            dateUpdated=date(2021, 11, 5),
        )
