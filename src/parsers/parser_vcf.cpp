#include <Variant.h>
#include "parsers/parser_vcf.hpp"

using namespace vcflib;
using namespace Anaquin;

bool ParserVCF::isVCF(const Reader &r)
{
    try
    {
        VariantCallFile vFile;
        
        auto src = std::string(r.src());
        vFile.open(src);
        
        if (!vFile.is_open())
        {
            return false;
        }
        
        vcflib::Variant var(vFile);

        for (auto i = 0; i < 5 && vFile.getNextVariant(var); i++) {}
        
        return true;
    }
    catch (...)
    {
        return false;
    }
}

