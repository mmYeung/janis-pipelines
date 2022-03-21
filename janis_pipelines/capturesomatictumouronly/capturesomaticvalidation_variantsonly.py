from janis_core import StringFormatter


from janis_pipelines.capturesomatictumouronly.capturesomatictumouronly_variantsonly import (
    CaptureSomaticTumourOnlyMultiCallersVariantsOnly,
)

from janis_bioinformatics.tools.variantcallers import (
    IlluminaSomaticPiscesVariantCallerTumourOnlyTargetedNoPON_5_2_10_49,
)


class CaptureSomaticValidationMultiCallersVariantsOnly(
    CaptureSomaticTumourOnlyMultiCallersVariantsOnly
):
    def id(self):
        return "CaptureSomaticValidationMultiCallersVariantsOnly"

    def friendly_name(self):
        return "Capture Somatic Tumour Only (Multi callers) [VARIANTS only]"

    def version(self):
        return "1.4.0"

    def add_pisces(self, bam_source):
        self.step(
            "vc_pisces",
            IlluminaSomaticPiscesVariantCallerTumourOnlyTargetedNoPON_5_2_10_49(
                bam=bam_source,
                sample_name=self.sample_name,
                referenceFolder=self.referenceFolder,
                intervals=self.intervals,
                minBQ=self.minBQ,
                minMQ=self.minMQ,
                minDP=self.minDP,
                minVAF=self.minVaf,
                VCminVQ=self.piscesVCminVQ,
                noise_level=self.pisces_noise_level,
                VQRminVQ=self.piscesVQRminVQ,
                pisces_awk_script=self.pisces_awk_script,
            ),
        )

        self.output(
            "pisces_bam",
            source=self.vc_pisces.out_bam,
            output_folder="Bam",
            output_name=StringFormatter(
                "{samplename}.pisces.bam", samplename=self.sample_name
            ),
        )
        self.output(
            "out_variants_pisces",
            source=self.vc_pisces.out,
            output_folder="variants",
            output_name=StringFormatter(
                "{samplename}.pisces.recode.vcf", samplename=self.sample_name
            ),
        )
        self.output(
            "variants_pisces",
            source=self.vc_pisces.variants,
            output_folder=["variants", "unfiltered"],
            output_name=StringFormatter(
                "{samplename}.pisces.vcf", samplename=self.sample_name
            ),
        )

