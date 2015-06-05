#include <string>
#include "resources.hpp"
#include "resources/silico.fa"
#include "resources/manual.txt"

/*
 * META Resources
 */

#include "resources/META.tab.fa"
#include "resources/META.ref.bed"
#include "resources/META.mix.csv"

/*
 * RNA Resources
 */

#include "resources/RNA.tab.fa"
#include "resources/RNA.ref.bed"
#include "resources/RNA.ref.gtf"
#include "resources/RNA.mix.csv"

/*
 * DNA Resources
 */

#include "resources/DNA.tab.fa"
#include "resources/DNA.ref.bed"
#include "resources/DNA.var.vcf"
#include "resources/DNA.mix.csv"

using namespace Spike;

ChromosomeInfo Resources::chromo()
{
    ChromosomeInfo i;
    
    i.v  = 1.0;
    i.id = "chrT";

    return i;
}

MixtureInfo Resources::mixture()
{
    MixtureInfo i;
    
    i.va = 1.0;
    i.vb = 1.0;

    return i;
}

#define ToString(x) std::string(reinterpret_cast<char*>(x))

std::string manual()
{
    return ToString(docs_manual_txt);
}

/*
 * Silico chromosome
 */

std::string silico_f()
{
    return ToString(data_silico_fa);
}

/*
 * META Resources
 */

std::string m_tab_f()
{
    return ToString(data_meta_META_tab_fa);
}

std::string m_mix_f()
{
    return ToString(data_meta_META_mix_csv);
}

std::string m_bed_f()
{
    return ToString(data_meta_META_ref_bed);
}

/*
 * RNA Resources
 */

std::string r_tab_f()
{
    return ToString(data_rna_RNA_tab_fa);
}

std::string r_gtf_f()
{
    return ToString(data_rna_RNA_ref_gtf);
}

std::string r_bed_f()
{
    return ToString(data_rna_RNA_ref_bed);
}

std::string r_mix_f()
{
    return ToString(data_rna_RNA_mix_csv);
}

/*
 * DNA Resources
 */

std::string d_tab_f()
{
    return ToString(data_dna_DNA_tab_fa);
}

std::string d_mix_f()
{
    return ToString(data_dna_DNA_mix_csv);
}

std::string d_bed_f()
{
    return ToString(data_dna_DNA_ref_bed);
}

std::string d_vcf_f()
{
    return ToString(data_dna_DNA_var_vcf);
}