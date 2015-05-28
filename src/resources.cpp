#include <string>
#include "resources.hpp"
#include "data/silico.fa"
#include "data/manual.txt"

/*
 * META Resources
 */

#include "data/META.mix.csv"

/*
 * RNA Resources
 */

#include "data/RNA.ref.bed"
#include "data/RNA.ref.gtf"
#include "data/RNA.mix.csv"

/*
 * DNA Resources
 */

#include "data/DNA.tab.fa"
#include "data/DNA.ref.bed"
#include "data/DNA.var.vcf"
#include "data/DNA.mix.csv"

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

std::string manual()
{
    return std::string(reinterpret_cast<char*>(docs_manual_txt));
}

/*
 * Silico chromosome
 */

std::string silico_f()
{
    return std::string(reinterpret_cast<char*>(data_silico_fa));
}

/*
 * META Resources
 */

std::string m_mix_f()
{
    return std::string(reinterpret_cast<char*>(data_meta_META_mix_csv));
}

/*
 * RNA Resources
 */

std::string r_gtf_f()
{
    return std::string(reinterpret_cast<char*>(data_rna_RNA_ref_gtf));
}

std::string r_bed_f()
{
    return std::string(reinterpret_cast<char*>(data_rna_RNA_ref_bed));
}

std::string r_mix_f()
{
    return std::string(reinterpret_cast<char*>(data_rna_RNA_mix_csv));
}

/*
 * DNA Resources
 */

std::string d_tab_f()
{
    return std::string(reinterpret_cast<char *>(data_dna_DNA_tab_fa));
}

std::string d_mix_f()
{
    return std::string(reinterpret_cast<char *>(data_dna_DNA_mix_csv));
}

std::string d_bed_f()
{
    return std::string(reinterpret_cast<char *>(data_dna_DNA_ref_bed));
}

std::string d_vcf_f()
{
    return std::string(reinterpret_cast<char *>(data_dna_DNA_var_vcf));
}