from janis_core import Int, Float, Directory, StringFormatter
from janis_bioinformatics.data_types import Bed, Vcf, VcfTabix


from janis_pipelines.capturesomatictumouronly_gatk.capturesomatictumouronlygatk_variantsonly import (
    CaptureSomaticTumourOnlyGATKVariantsOnly,
)

from janis_bioinformatics.tools.variantcallers import (
    VardictGermlineVariantCaller,
    VarscanGermlineCNSVariantCaller,
    IlluminaSomaticPiscesVariantCallerTumourOnlyTargeted_5_2_10_49,
)

from janis_bioinformatics.tools.dawson import (
    GenerateChromosomeIntervalsFromBed,
    GenerateVarscan2HeaderLines,
)

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
        self.input("referenceFolder", Directory())

    def add_inputs_for_intervals(self):
        super().add_inputs_for_intervals()
        self.input("intervals", Bed())

    def add_inputs_for_vc(self):
        # General parameters
        self.input("minVaf", Float())
        self.input("minMQ", Int())
        self.input("minAD", Int())
        self.input("minDP", Int())
        self.input("minBQ", Int())
        # Varscan2
        self.input("pileupMaxDepth", Int(), default=20000)
        self.input("pileupMinBQ", Int())
        self.input("varscanPval", Float(), default=0.0001)

        # Pisces
        self.input("piscesVCminVQ", Int())
        self.input("piscesVQRminVQ", Int())

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
                allele_freq_threshold=self.minVaf,
                minMappingQual=self.minMQ,
                sv=False,
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
                maxDepth=self.pileupMaxDepth,
                minBQ=self.minBQ,
                minDP=self.minDP,
                minAD=self.minAD,
                minVaf=self.minVaf,
                pval=self.varscanPval,
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
                referenceFolder=self.referenceFolder,
                PON=self.panel_of_normals,
                intervals=self.intervals,
                minBQ=self.minBQ,
                VCminVQ=self.piscesVCminVQ,
                VQRminVQ=self.piscesVQRminVQ,
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

    def add_inputs_for_intervals(self):
        self.input("intervals", Bed())

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

