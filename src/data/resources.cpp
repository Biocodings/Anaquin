#include <string>

/*
 * Scripts
 */

#include "resources/linear.R"
#include "resources/manual.txt"

/*
 * Fusion Resources
 */

#include "resources/FUS.v1.ref"
#include "resources/FUS.v3.csv"
#include "resources/FusionGenes.chrTv1.bed"
#include "resources/NormalParentGenes.chrTv1.bed"

/*
 * Conjoint Resources
 */

#include "resources/Ladder_v3.csv"

/*
 * RNA Resources
 */

#include "resources/RNA_1.gtf"
#include "resources/RNA_4_1.csv"

/*
 * META Resources
 */

#include "resources/META_v6.csv"
#include "resources/META_v1_tab.fa"
#include "resources/META_v1_tab.bed"

/*
 * DNA Resources
 */

#include "resources/DNA_v3.csv"
#include "resources/DNA.variant.bed"
#include "resources/DNA.standards.chrT.gtf"

#define ToString(x) std::string(reinterpret_cast<char*>(x))

std::string ChromoName()
{
    return "chrT";
}

float ChromoVersion()
{
    return 1.0;
}

float MixtureVersion()
{
    return 1.0;
}

std::string LinearR()
{
    return ToString(data_linear_R);
}

std::string Manual()
{
    return ToString(data_manual_txt);
}

/*
 * Fusion Resources
 */

std::string FusionDataMix()
{
    return ToString(data_fusion_FUS_v3_csv);
}

std::string FusionDataRef()
{
    return ToString(data_fusion_FUS_v1_ref);
}

std::string FusionMutatedRef()
{
    return ToString(data_fusion_FusionGenes_chrTv1_bed);
}

std::string FusionNormalRef()
{
    return ToString(data_fusion_NormalParentGenes_chrTv1_bed);
}

/*
 * Ladder Resources
 */

std::string LadderDataMix()
{
    return ToString(data_ladder_Ladder_v3_csv);
}

/*
 * META Resources
 */

std::string MetaDataMix()
{
    return ToString(data_meta_META_v6_csv);
}

std::string MetaDataBed()
{
    return ToString(data_meta_META_v1_tab_bed);
}

/*
 * Transcriptome Resources
 */

std::string TransStandGTF()
{
    return ToString(data_trans_RNA_1_gtf);
}

std::string TransDataMix()
{
    return ToString(data_trans_RNA_4_1_csv);
}

/*
 * Variant Resources
 */

std::string VarDataMix()
{
    return ToString(data_var_DNA_v3_csv);
}

std::string VarDataBed()
{
    return ToString(data_var_DNA_variant_bed);
}

std::string VarStandGTF()
{
    return ToString(data_var_DNA_standards_chrT_gtf);
}