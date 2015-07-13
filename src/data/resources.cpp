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

/*
 * Conjoint Resources
 */

#include "resources/CON.v3.mix.csv"

/*
 * RNA Resources
 */

//#include "resources/chrT.v2.fa"
#include "resources/RNA.v1.bed"
#include "resources/RNA.v1.gtf"
#include "resources/RNA.v4.1.mix"
#include "resources/RNA.usage.txt"

/*
 * META Resources
 */

#include "resources/META.usage.txt"
#include "resources/META.v1.tab.fa"
#include "resources/META.v6.mix.csv"
#include "resources/META.v1.tab.bed"

/*
 * DNA Resources
 */

#include "resources/DNA.usage.txt"
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

std::string FusDataMix()
{
    return ToString(data_fus_FUS_v3_csv);
}

std::string FusDataRef()
{
    return ToString(data_fus_FUS_v1_ref);
}

/*
 * Ladder Resources
 */

std::string LadderDataMix()
{
    return ToString(data_con_CON_v3_mix_csv);
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
 * RNA Resources
 */

//std::string RNADataFA()
//{
  //  return ToString(data_rna_chrT_v2_fa);
//}

std::string RNADataGTF()
{
    return ToString(data_rna_RNA_v1_gtf);
}

std::string RNADataBed()
{
    return ToString(data_rna_RNA_v1_bed);
}

std::string RNADataMix()
{
    return ToString(data_rna_RNA_v4_1_mix);
}

/*
 * Variant Resources
 */

std::string VARDataMix()
{
    return ToString(data_dna_DNA_v3_mix_csv);
}

std::string VARDataBed()
{
    return ToString(data_dna_DNA_variant_bed);
}