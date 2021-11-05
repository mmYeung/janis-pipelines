from janis_core import String, Array
from janis_unix.tools import UncompressArchive
from janis_core.operators.standard import FirstOperator
from janis_core import WorkflowMetadata

from janis_bioinformatics.data_types import BamBai, Vcf, CompressedVcf

from janis_pipelines.wgs_somatic_gatk.wgssomaticgatk_variantsonly import (
    WGSSomaticGATKVariantsOnly,
)
from janis_bioinformatics.tools.pmac import GenerateIntervalsByChromosome
from janis_bioinformatics.tools.common.gatkbasecalbam import (
    GATKBaseRecalBQSRWorkflow_4_1_3,
)
from janis_bioinformatics.tools.variantcallers import (
    GatkSomaticVariantCallerTumorOnlyTargeted,
)
from janis_bioinformatics.tools.bcftools import (
    BcfToolsConcat_1_9,
    BcfToolsSort_1_9,
)
from janis_bioinformatics.tools.pmac import AddBamStatsGermline_0_1_0


class CaptureSomaticTumourOnlyGATKVariantsOnly(WGSSomaticGATKVariantsOnly):
    def id(self):
        return "CaptureSomaticGATKVariantsOnly"

    def friendly_name(self):
        return "Capture Somatic Tumour Only(GATK only) [VARIANTS only]"

    def version(self):
        return "1.4.0"

    def add_inputs(self):
        self.input("tumour_bam", BamBai())
        self.input("tumour_name", String())
        self.add_inputs_for_configuration()
        self.add_inputs_for_intervals()
        self.add_inputs_for_reference()

    def add_gatk_variantcaller(self, tumour_bam_source):
        if "generate_gatk_intervals" in self.step_nodes:
            generated_intervals = self.generate_gatk_intervals.out_regions
        else:
            generated_intervals = self.step(
                "generate_gatk_intervals",
                GenerateIntervalsByChromosome(reference=self.reference),
                when=self.gatk_intervals.is_null(),
            ).out_regions

        intervals = FirstOperator([self.gatk_intervals, generated_intervals])

        recal_ins = {
            "reference": self.reference,
            "intervals": intervals,
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
                intervals=intervals,
                reference=self.reference,
                gnomad=self.gnomad,
                panel_of_normals=self.panel_of_normals,
            ),
            scatter=["intervals", "bam"],
        )

        self.step(
            "vc_gatk_merge",
            BcfToolsConcat_1_9(vcf=self.vc_gatk.out.as_type(Array(Vcf))),
        )
        self.step(
            "vc_gatk_sort_combined",
            BcfToolsSort_1_9(vcf=self.vc_gatk_merge.out.as_type(CompressedVcf)),
        )

        self.step(
            "vc_gatk_uncompressvcf",
            UncompressArchive(file=self.vc_gatk_sort_combined.out),
        )

    def add_addbamstats(self, tumour_bam_source):
        self.step(
            "addbamstats",
            AddBamStatsGermline_0_1_0(
                bam=tumour_bam_source,
                vcf=self.vc_gatk_uncompressvcf.out.as_type(Vcf),
                reference=self.reference,
            ),
        )

    def constructor(self):
        self.add_inputs()

        self.add_gatk_variantcaller(tumour_bam_source=self.tumour_bam)

        self.add_addbamstats(tumour_bam_source=self.tumour_bam)

        self.output(
            "out_variants_gatk",
            source=self.vc_gatk_sort_combined.out,
            output_folder="variants",
            doc="Merged variants from the GATK caller",
        )
        self.output(
            "out_variants_gakt_split",
            source=self.vc_gatk.out,
            output_folder=["variants", "byInterval"],
            doc="Unmerged variants from the GATK caller (by interval)",
        )
        self.output(
            "out_variants_bamstats",
            source=self.addbamstats.out,
            output_folder="variants",
            doc="Final vcf",
        )

    def bind_metadata(self):
        from datetime import date

        return WorkflowMetadata(
            version="1.0.0",
            contributors=["Miriam M Yeung"],
            dateCreated=date(2021, 11, 5),
            dateUpdated=date(2021, 11, 5),
        )
