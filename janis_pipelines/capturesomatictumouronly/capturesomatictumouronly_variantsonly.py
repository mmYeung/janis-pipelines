from janis_core import Int, Float, String, Directory, File, StringFormatter
from janis_bioinformatics.data_types import Vcf


from janis_pipelines.capturesomatictumouronly_gatk.capturesomatictumouronlygatk_variantsonly import (
    CaptureSomaticTumourOnlyGATKVariantsOnly,
)

from janis_bioinformatics.tools.variantcallers import (
    VardictGermlineVariantCaller,
    VarscanGermlineCNSVariantCaller,
    IlluminaSomaticPiscesVariantCallerTumourOnlyTargeted_5_2_10_49,
)

from janis_bioinformatics.tools.dawson import GenerateVarscan2HeaderLines

from janis_bioinformatics.tools.pmac import (
    CombineVariants_0_0_8,
    AddBamStatsGermline_0_1_0,
    GenerateVardictHeaderLines,
)

from janis_bioinformatics.tools.bcftools import BcfToolsSort_1_9

from janis_unix.tools import UncompressArchive

from janis_bioinformatics.tools.htslib import BGZipLatest, TabixLatest


class CaptureSomaticTumourOnlyMultiCallersVariantsOnly(
    CaptureSomaticTumourOnlyGATKVariantsOnly
):
    def id(self):
        return "CaptureSomaticTumourOnlyMultiCallersVariantsOnly"

    def friendly_name(self):
        return "Capture Somatic Tumour Only (Multi callers) [VARIANTS only]"

    def version(self):
        return "1.4.0"

    def constructor(self):
        self.add_inputs()

        self.add_gatk_variantcaller(tumour_bam_source=self.bam)
        ## Adding outputs from GATK for easier referencing in following pipelines
        self.output(
            "out_variants_gatk",
            source=self.vc_gatk_filterpass.out,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}.gatk.recode", samplename=self.sample_name
            ),
        )
        self.output(
            "variants_gatk",
            source=self.vc_gatk_sort.out,
            output_folder=["variants", "unfiltered"],
            output_name=StringFormatter(
                "{samplename}.gatk", samplename=self.sample_name
            ),
        )
        self.output(
            "gatk_bam",
            source=self.vc_gatk.out_bam,
            output_folder=["Bam"],
            output_name=StringFormatter(
                "{samplename}.gatk", samplename=self.sample_name
            ),
        )

        self.add_vardict(bam_source=self.bam)
        self.add_varscan2(bam_source=self.bam)
        self.add_pisces(bam_source=self.bam)

        self.add_combine_variants(bam_source=self.bam)

    def add_inputs_for_reference(self):
        super().add_inputs_for_reference()
        self.input("reference_folder", Directory())

    def add_inputs_for_vc(self):
        # General parameters
        self.input("min_vaf", Float())
        self.input("min_mq", Int())
        self.input("min_ad", Int())
        self.input("min_dp", Int())
        self.input("min_bq", Int())
        # Varscan2
        self.input("pileup_max_depth", Int(), default=20000)
        self.input("pileup_min_bq", Int())
        self.input("varscan_pval", Float(), default=0.0001)

        # Pisces
        self.input("pisces_ploidy", String(optional=True))
        self.input("pisces_vc_min_vq", Int())
        self.input("pisces_vqr_min_vq", Int())
        self.input("pisces_noise_level", Int(optional=True))
        self.input("pisces_awk_script", File())

    def add_inputs(self):
        super().add_inputs()
        self.add_inputs_for_vc()

    def add_vardict(self, bam_source):
        self.step(
            "generate_vardict_headerlines",
            GenerateVardictHeaderLines(reference=self.reference),
        )

        self.step(
            "vc_vardict",
            VardictGermlineVariantCaller(
                bam=bam_source,
                sample_name=self.sample_name,
                intervals=self.intervals,
                header_lines=self.generate_vardict_headerlines.out,
                reference=self.reference,
                allele_freq_threshold=self.min_vaf,
                min_mapping_qual=self.min_mq,
                no_sv_call=True,
            ),
        )

        self.output(
            "out_variants_vardict",
            source=self.vc_vardict.out,
            output_folder="variants",
            output_name=StringFormatter(
                "{samplename}.vardict.recode.vcf", samplename=self.sample_name
            ),
        )
        self.output(
            "variants_vardict",
            source=self.vc_vardict.variants,
            output_folder=["variants", "unfiltered"],
            output_name=StringFormatter(
                "{samplename}.vardict.vcf", samplename=self.sample_name
            ),
        )

    def add_varscan2(self, bam_source):
        self.step(
            "generate_varscan2_headerlines",
            GenerateVarscan2HeaderLines(reference=self.reference),
        )
        self.step(
            "vc_varscan2",
            VarscanGermlineCNSVariantCaller(
                sample_name=self.sample_name,
                bam=bam_source,
                reference=self.reference,
                pileup_max_depth=self.pileup_max_depth,
                min_bq=self.min_bq,
                min_dp=self.min_dp,
                min_ad=self.min_ad,
                min_vaf=self.min_vaf,
                pval=self.varscan_pval,
                header_lines=self.generate_varscan2_headerlines.out,
            ),
        )

        self.output(
            "out_variants_varscan2",
            source=self.vc_varscan2.out,
            output_folder="variants",
            output_name=StringFormatter(
                "{samplename}.varscan.recode.vcf", samplename=self.sample_name
            ),
        )
        self.output(
            "variants_varscan2",
            source=self.vc_varscan2.variants,
            output_folder=["variants", "unfiltered"],
            output_name=StringFormatter(
                "{samplename}.varscan.vcf", samplename=self.sample_name
            ),
        )

    def add_pisces(self, bam_source):
        self.step(
            "vc_pisces",
            IlluminaSomaticPiscesVariantCallerTumourOnlyTargeted_5_2_10_49(
                bam=bam_source,
                sample_name=self.sample_name,
                reference_folder=self.reference_folder,
                pon=self.panel_of_normals,
                intervals=self.intervals,
                ploidy=self.pisces_ploidy,
                min_bq=self.min_bq,
                min_mq=self.min_mq,
                min_vaf=self.min_vaf,
                vc_min_vq=self.pisces_vc_min_vq,
                vqr_min_vq=self.pisces_vqr_min_vq,
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

    def add_combine_variants(self, bam_source):
        self.step(
            "combine_variants",
            CombineVariants_0_0_8(
                vcfs=[
                    self.vc_gatk_filterpass.out,
                    self.vc_pisces.out,
                    self.vc_varscan2.out,
                    self.vc_vardict.out,
                ],
                type="germline",
                columns=["AC", "AN", "AF", "AD", "DP", "GT"],
            ),
        )

        self.step(
            "combined_sort", BcfToolsSort_1_9(vcf=self.combine_variants.out)
        )

        self.step(
            "combined_uncompress",
            UncompressArchive(file=self.combined_sort.out),
        )

        self.step(
            "combined_addbamstats",
            AddBamStatsGermline_0_1_0(
                bam=bam_source,
                vcf=self.combined_uncompress.out.as_type(Vcf),
                reference=self.reference,
            ),
        )

        self.step(
            "compress_combined_addbamstats",
            BGZipLatest(file=self.combined_addbamstats),
        )

        self.step(
            "tabixvcf", TabixLatest(inp=self.compress_combined_addbamstats.out)
        )

        self.output(
            "out_variants",
            source=self.tabixvcf.out,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{samplename}.combined.vcf", samplename=self.sample_name
            ),
            doc="Combined variants from all 4 variant callers",
        )

