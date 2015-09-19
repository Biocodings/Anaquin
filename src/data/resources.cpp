#include <string>

/*
 * Scripts
 */

#include "resources/linear.R"
#include "resources/viewer.py"
#include "resources/manual.txt"

/*
 * Fusion Resources
 */

#include "resources/FUSBreak_1.0.ref"
#include "resources/FUSMixture_3.0.csv"
#include "resources/FUSFusionStandard_1.0.bed"
#include "resources/FUSNormalStandard_1.0.bed"

/*
 * Conjoint Resources
 */

#include "resources/LadderMixture_3.0.csv"

/*
 * RNA Resources
 */

#include "resources/MTR002.v013.csv"
#include "resources/MTR003.v013.csv"
#include "resources/MTR004.v013.csv"
#include "resources/ATR001.v032.gtf"

/*
 * Metagenomics Resources
 */

#include "resources/META_v1_tab.fa"
#include "resources/MME023.v013.csv"
#include "resources/METAStandard_1.0.bed"

/*
 * Variant Resources
 */

#include "resources/AVA009.v032.vcf"
#include "resources/MVA011.v013.csv"
#include "resources/MVA012.v013.csv"
#include "resources/AVA008.v032.bed"

#define ToString(x) std::string(reinterpret_cast<char*>(x))

std::string ViewerScript()
{
    return ToString(scripts_viewer_py);
}

std::string Manual()
{
    return ToString(data_manual_txt);
}

/*
 * R script for plotting
 */

std::string RScriptCoverage()
{
    return ToString(data_linear_R);
}

/*
 * Fusion Resources
 */

std::string FusionDataMix()
{
    return ToString(data_fusion_FUSMixture_3_0_csv);
}

std::string FusionDataRef()
{
    return ToString(data_fusion_FUSBreak_1_0_ref);
}

/*
std::string FusionMutatedRef()
{
    return ToString(data_fusion_FusionGenes_chrTv1_bed);
}

std::string FusionNormalRef()
{
    return ToString(data_fusion_NormalParentGenes_chrTv1_bed);
}
*/

/*
 * Ladder Resources
 */

std::string LadderDataMix()
{
    return ToString(data_ladder_LadderMixture_3_0_csv);
}

/*
 * META Resources
 */

std::string MetaDataMix()
{
    return ToString(data_meta_MME023_v013_csv);
}

std::string MetaDataBed()
{
    return ToString(data_meta_METAStandard_1_0_bed);
}

/*
 * Transcriptome Resources
 */

std::string TransStandGTF()
{
    return ToString(data_trans_ATR001_v032_gtf);
}

std::string TransDataMixA()
{
    return ToString(data_trans_MTR002_v013_csv);
}

std::string TransDataMixB()
{
    return ToString(data_trans_MTR003_v013_csv);
}

std::string TransDataMixAB()
{
    return ToString(data_trans_MTR004_v013_csv);
}

/*
 * Variant Resources
 */

std::string VarDataMixA()
{
    return ToString(data_var_MVA011_v013_csv);
}

std::string VarDataMixF()
{
    return ToString(data_var_MVA012_v013_csv);
}

std::string VarDataVCF()
{
    return ToString(data_var_AVA009_v032_vcf);
}

std::string VarDataBed()
{
    return ToString(data_var_AVA008_v032_bed);
}