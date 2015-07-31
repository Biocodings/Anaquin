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
#include "resources/fusion.bed"
#include "resources/normal.bed"

/*
 * Conjoint Resources
 */

#include "resources/CON.v3.mix.csv"

/*
 * RNA Resources
 */

#include "resources/RNA.v1.gtf"
#include "resources/RNA.v4.1.csv"

/*
 * META Resources
 */

#include "resources/META.v1.tab.fa"
#include "resources/META.v6.mix.csv"
#include "resources/META.v1.tab.bed"

/*
 * DNA Resources
 */

#include "resources/DNA.v3.mix.csv"
#include "resources/DNA.variant.bed"

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
    return ToString(data_fusion_fusion_bed);
}

std::string FusionNormalRef()
{
    return ToString(data_fusion_normal_bed);
}

/*
 * Ladder Resources
 */

std::string LadderDataMix()
{
    return ToString(data_ladder_CON_v3_mix_csv);
}

/*
 * META Resources
 */

std::string MetaDataTab()
{
    return ToString(data_meta_META_v1_tab_fa);
}

std::string MetaDataMix()
{
    return ToString(data_meta_META_v6_mix_csv);
}

std::string MetaDataBed()
{
    return ToString(data_meta_META_v1_tab_bed);
}

/*
 * Transcriptome Resources
 */

std::string TransDataGTF()
{
    return ToString(data_trans_RNA_v1_gtf);
}

std::string TransDataMix()
{
    return ToString(data_trans_RNA_v4_1_csv);
}

/*
 * Variant Resources
 */

std::string VarDataMix()
{
    return ToString(data_var_DNA_v3_mix_csv);
}

std::string VarDataBed()
{
    return ToString(data_var_DNA_variant_bed);
}