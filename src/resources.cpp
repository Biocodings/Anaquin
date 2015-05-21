#include <string>
#include "data/silico.fa"
#include "data/manual.txt"
#include "data/standards.bed"
#include "data/standards.gtf"
#include "data/rna_standards.txt"

using namespace std;

std::string silico_fa()
{
    return std::string(reinterpret_cast<char*>(data_silico_fa));
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