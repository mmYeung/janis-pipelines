from janis_core import String, Array, WorkflowMetadata, StringFormatter
from janis_unix.tools import UncompressArchive

from janis_bioinformatics.data_types import BamBai, Bed, Vcf, CompressedVcf

from janis_pipelines.wgs_somatic_gatk.wgssomaticgatk_variantsonly import (
    WGSSomaticGATKVariantsOnly,
)

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


class CaptureSomaticTumourOnlyGATKVariantsOnly(WGSSomaticGATKVariantsOnly):
    def id(self):
        return "CaptureSomaticGATKVariantsOnly"

    def friendly_name(self):
        return "Capture Somatic Tumour Only(GATK only) [VARIANTS only]"

    def version(self):
        return "1.4.0"

    def add_inputs(self):
        self.input("bam", BamBai())
        self.input("sample_name", String())
        self.input("intervals", Bed())
        self.add_inputs_for_configuration()
        self.add_inputs_for_intervals()
        self.add_inputs_for_reference()

    def add_gatk_variantcaller(self, tumour_bam_source):

        # intervals = FirstOperator([self.gatk_intervals, generated_intervals])

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
            scatter="intervals",
        )

        self.step(
            "vc_gatk",
            GatkSomaticVariantCallerTumorOnlyTargeted(
                bam=self.bqsr_tumour.out,
                intervals=self.intervals,
                reference=self.reference,
                gnomad=self.gnomad,
                panel_of_normals=self.panel_of_normals,
            ),
            scatter=["intervals", "bam"],
        )

        self.step(
            "vc_gatk_merge",
            BcfToolsConcat_1_9(vcf=self.vc_gatk.variants.as_type(Array(Vcf))),
        )
        self.step(
            "vc_gatk_sort_combined",
            BcfToolsSort_1_9(vcf=self.vc_gatk_merge.out.as_type(CompressedVcf)),
        )

        self.step(
            "vc_gatk_normalise",
            BcfToolsNorm_1_9(vcf=self.vc_gatk_sort_combined.out),
        )

        self.step(
            "vc_gatk_uncompressvcf",
            UncompressArchive(file=self.vc_gatk_sort_combined.out),
        )

        self.step(
            "filterpass",
            VcfToolsvcftools_0_1_16(
                vcf=self.uncompress.out.as_type(Vcf),
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
            source=self.vc_gatk_sort_combined.out,
            output_folder=["variants", "unfiltered"],
            doc="Merged variants from the GATK caller",
        )

        self.output(
            "out_variants_pass_gatk",
            source=self.filterpass.out,
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
