//#include "parsers/parser_vcf.hpp"
//#include <Variant.h>
//
//using namespace Anaquin;
//
//bool ParserVCF::isVCF(const Reader &r)
//{
//    try
//    {
//        vcflib::VariantCallFile vFile;
//        
//        auto src = std::string(r.src());
//        vFile.open(src);
//        
//        if (!vFile.is_open())
//        {
//            return false;
//        }
//        
//        vcflib::Variant var(vFile);
//
//        for (auto i = 0; i < 5 && vFile.getNextVariant(var); i++) {}
//        
//        return true;
//    }
//    catch (...)
//    {
//        return false;
//    }
//}
//
