from janis_core import (
    Int,
    Float,
    String,
    Array,
    Directory,
    File,
    StringFormatter,
    ScatterMethod,
    ScatterDescription,
)

from janis_bioinformatics.tools import BioinformaticsWorkflow

from janis_bioinformatics.data_types import BamBai, Bed, VcfTabix, FastaWithDict

from janis_pipelines.capturesomatictumouronly.capturesomatictumouronly_variantsonly import (
    CaptureSomaticTumourOnlyMultiCallersVariantsOnly,
)


class CaptureSomaticsTumourOnlyVariantsOnlyMultiSample(BioinformaticsWorkflow):
    def id(self):
        return "CaptureSomaticsTumourOnlyVariantsOnlyMultiSample"

    def friendly_name(self):
        return "Capture Somatic Tumour Only Variants Only (MultiSample)"

    def version(self):
        return "0.0.0"

    def constructor(self):
        ##Inputs
        self.input("bams", Array(BamBai()))
        self.input("sample_name", Array(String()))
        ## Intervals
        self.input("intervals", Bed())
        # References
        self.input("reference", FastaWithDict())
        self.input("reference_alt", File())
        self.input("reference_folder", Directory())
        # GATK references
        self.input("snps_dbsnp", VcfTabix())
        self.input("snps_1000gp", VcfTabix())
        self.input("known_indels", VcfTabix())
        self.input("mills_indels", VcfTabix())
        self.input("gnomad", VcfTabix())
        self.input("panel_of_normals", VcfTabix())
        # Variant Callers General
        self.input("min_bq", Int())
        self.input("min_mq", Int())
        self.input("min_ad", Int())
        self.input("min_dp", Int())
        self.input("min_vaf", Float())
        ## Varscan2
        self.input("pileup_min_bq", Int())
        self.input("pileup_max_depth", Int(), default=20000)
        self.input("varscan_pval", Float(), default=0.0001)
        ## Pisces
        self.input("pisces_ploidy", String(optional=True))
        self.input("pisces_vc_min_vq", Int())
        self.input("pisces_vqr_min_vq", Int())
        self.input("pisces_noise_level", Int(optional=True))
        self.input("pisces_awk_script", File())

        self.step(
            "capture_variantcalling",
            CaptureSomaticTumourOnlyMultiCallersVariantsOnly(
                bam=self.bams,
                sample_name=self.sample_name,
                intervals=self.intervals,
                reference=self.reference,
                reference_folder=self.reference_folder,
                known_indels=self.known_indels,
                mills_indels=self.mills_indels,
                snps_1000gp=self.snps_1000gp,
                snps_dbsnp=self.snps_dbsnp,
                gnomad=self.gnomad,
                panel_of_normals=self.panel_of_normals,
                min_bq=self.min_bq,
                min_mq=self.min_mq,
                min_ad=self.min_ad,
                min_dp=self.min_dp,
                min_vaf=self.min_vaf,
                pileup_min_bq=self.pileup_min_bq,
                pileup_max_depth=self.pileup_max_depth,
                varscan_pval=self.varscan_pval,
                pisces_ploidy=self.pisces_ploidy,
                pisces_vc_min_vq=self.pisces_vc_min_vq,
                pisces_vqr_min_vq=self.pisces_vqr_min_vq,
                pisces_noise_level=self.pisces_noise_level,
                pisces_awk_script=self.pisces_awk_script,
            ),
            scatter=ScatterDescription(
                ["bam", "sample_name"],
                method=ScatterMethod.dot,
                labels=self.sample_name,
            ),
        )
        self.capture_outputs_from_step(self.capture_variantcalling)
