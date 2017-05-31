#include "tools/tools.hpp"
#include "VarQuin/v_structure.hpp"

using namespace Anaquin;

typedef SeqVariant::Context Context;

extern Scripts PlotVLODR();
extern Scripts PlotVGROC();
extern Scripts PlotVCROC();
extern Scripts PlotAllele();

extern Path __output__;
extern std::string __full_command__;

VStructure::EStats VStructure::analyzeE(const FileName &file, const Options &o)
{
    const auto regs = Standard::instance().r_var.regs1();
    
    VStructure::EStats stats;
    
    if (!file.empty())
    {
        readVFile(file, [&](const ParserVCF::Data &x, const ParserProgress &p)
        {
            // Only interested whether the variant falls in the reference regions
            if (regs.count(x.cID) && regs.at(x.cID).contains(x.l))
            {
                stats.found++;
            }
        });
    }
    else
    {
        stats.found = NAN;
    }

    return stats;
}

VStructure::SStats VStructure::analyzeS(const FileName &file, const Options &o)
{
    VStructure::SStats stats;
    return stats;
}

//static void writeQuins(const FileName &file,
//                       const VStructure::SStats &ss,
//                       const VStructure::Options &o)
//{
//}
//
//static void writeDetected(const FileName &file,
//                          const VStructure::SStats &ss,
//                          const VStructure::Options &o)
//{
//  
//}
//
//static void writeSummary(const FileName &file,
//                         const FileName &endo,
//                         const FileName &seqs,
//                         const VStructure::EStats &es,
//                         const VStructure::SStats &ss,
//                         const VStructure::Options &o)
//{
//   
//}

void VStructure::report(const FileName &endo, const FileName &seqs, const Options &o)
{
   
}
