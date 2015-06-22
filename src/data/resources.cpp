#include <string>

/*
 * Scripts
 */

#include "resources/linear.R"
#include "resources/manual.txt"

/*
 * RNA Resources
 */

#include "resources/r_manual.txt"
#include "resources/chrT.v1.fa"
#include "resources/RNA.v1.bed"
#include "resources/RNA.v1.gtf"
#include "resources/RNA.v4.1.mix"
#include "resources/RNA.v1.2.tab.fa"

/*
 * META Resources
 */

#include "resources/m_manual.txt"
#include "resources/META.v1.tab.fa"
#include "resources/META.v6.mix.csv"
#include "resources/META.v1.tab.bed"

/*
 * DNA Resources
 */

#include "resources/d_manual.txt"
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

std::string RNADataFA()
{
    return ToString(data_chrT_v1_fa);
}

std::string RNADataTab()
{
    return ToString(data_rna_RNA_v1_2_tab_fa);
}

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
 * DNA Resources
 */

std::string DNADataMix()
{
    return ToString(data_dna_DNA_v3_mix_csv);
}

std::string DNADataBed()
{
    return ToString(data_dna_DNA_variant_bed);
}