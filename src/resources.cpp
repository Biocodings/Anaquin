#include <string>
#include "data/silico.fa"
#include "data/manual.txt"
#include "data/DNA.tab.fa"
#include "data/DNA.ref.bed"
#include "data/standards.bed"
#include "data/standards.gtf"
#include "data/rna_standards.txt"
#include "data/variant.ChrT51.vcf"

using namespace std;

std::string silico_fa()
{
    return std::string(reinterpret_cast<char*>(data_silico_fa));
}

std::string r_standards_txt()
{
    return std::string(reinterpret_cast<char*>(data_rna_rna_standards_txt));
}

std::string d_ref_bed()
{
    return std::string(reinterpret_cast<char*>(data_dna_DNA_ref_bed));
}

std::string d_tab_fa()
{
    return std::string(reinterpret_cast<char*>(data_dna_DNA_tab_fa));
}

std::string d_variant_vcf()
{
    return std::string(reinterpret_cast<char*>(data_dna_variant_ChrT51_vcf));
}

std::string standards_gtf()
{
    return std::string(reinterpret_cast<char*>(data_rna_standards_gtf));
}

std::string standards_bed()
{
    return std::string(reinterpret_cast<char*>(data_rna_standards_bed));
}

std::string rna_standards_txt()
{
    return std::string(reinterpret_cast<char*>(data_rna_rna_standards_txt));
}

std::string manual()
{
    return std::string(reinterpret_cast<char*>(docs_manual_txt));
}