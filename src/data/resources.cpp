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

#include "resources/AFU004.v032.ref"
#include "resources/AFU005.v032.bed"
#include "resources/MFU007.v013.csv"

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

#include "resources/MME023.v013.csv"
#include "resources/AME015.v032.bed"

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

std::string FusionDataMixA()
{
    return ToString(data_fusion_MFU007_v013_csv);
}

std::string FusionDataRef()
{
    return ToString(data_fusion_AFU004_v032_ref);
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
    return ToString(data_meta_AME015_v032_bed);
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