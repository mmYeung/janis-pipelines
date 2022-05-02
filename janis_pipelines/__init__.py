from .__meta__ import __version__

# from janis_pipelines.alignment.alignment import BwaAlignment

from janis_pipelines.wgs_germline.wgsgermline import (
    WGSGermlineMultiCallers,
    WGSGermlineMultiCallersVariantsOnly,
)
from janis_pipelines.wgs_germline_gatk.wgsgermlinegatk import (
    WGSGermlineGATK,
    WGSGermlineGATKVariantsOnly,
)

from janis_pipelines.wgs_somatic.wgssomatic import (
    WGSSomaticMultiCallers,
    WGSSomaticMultiCallersVariantsOnly,
)
from janis_pipelines.wgs_somatic_gatk.wgssomaticgatk import (
    WGSSomaticGATK,
    WGSSomaticGATKVariantsOnly,
)

from janis_pipelines.capturesomatictumouronly_gatk.capturesomatictumouronlygatk_variantsonly import (
    CaptureSomaticTumourOnlyGATKVariantsOnly,
)

from janis_pipelines.capturesomatictumouronly_gatk.capturesomaticvalidationgatk_variantsonly import (
    CaptureSomaticValidationGATKVariantsOnly,
)

from janis_pipelines.capturesomatictumouronly.capturesomatictumouronly import (
    CaptureSomaticTumourOnly,
)

from janis_pipelines.capturesomatictumouronly.capturesomatictumouronly_variantsonly import (
    CaptureSomaticTumourOnlyMultiCallersVariantsOnly,
)

from janis_pipelines.capturesomatictumouronly.capturesomatictumouronly_variantsonlymultisample import (
    CaptureSomaticsTumourOnlyVariantsOnlyMultiSample,
)

from janis_pipelines.capturesomatictumouronly.capturesomaticvalidation_variantsonly import (
    CaptureSomaticValidationMultiCallersVariantsOnly,
)

from janis_pipelines.capturesomatictumouronly.capturesomaticvalidation_variantsonlymultisample import (
    CaptureSomaticValidationVariantsOnlyMultiSample,
)

from janis_pipelines.capturesomatictumouronly.capturesomatictumouronlyumi import (
    CaptureSomaticTumourOnlyUMI,
)

from janis_pipelines.capturesomatictumouronly.capturesomaticvalidationumi import (
    CaptureSomaticValidationUMI,
)

from janis_pipelines.capturesomatictumouronly.capturesomatictumouronlyumimulitsample import (
    CaptureSomaticTumourOnlyUMIMultiSample,
)


from janis_pipelines.capturesomatictumouronly.capturesomaticvalidationumimulitsample import (
    CaptureSomaticValidationUMIMultiSample,
)

from janis_pipelines.alignment.agilentUMIalignment import AgilentUMIalignment
