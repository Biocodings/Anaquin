#include "tools/tools.hpp"
#include "VarQuin/v_conjoint.hpp"

using namespace Anaquin;

typedef SequinVariant::Context Context;

extern Scripts PlotVLODR();
extern Scripts PlotVGROC();
extern Scripts PlotVCROC();
extern Scripts PlotAllele();

extern Path __output__;
extern std::string __full_command__;

VConjoint::EStats VConjoint::analyzeE(const FileName &file, const Options &o)
{
    const auto regs = Standard::instance().r_var.regs1();
    
    VConjoint::EStats stats;
    
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

VConjoint::SStats VConjoint::analyzeS(const FileName &file, const Options &o)
{
//    const auto &r = Standard::instance().r_var;

    VConjoint::SStats stats;
    
    return stats;
}

//static void writeQuins(const FileName &file,
//                       const VConjoint::SStats &ss,
//                       const VConjoint::Options &o)
//{
//
//}
//
//static void writeDetected(const FileName &file,
//                          const VConjoint::SStats &ss,
//                          const VConjoint::Options &o)
//{
//}
//
//static void writeSummary(const FileName &file,
//                         const FileName &endo,
//                         const FileName &seqs,
//                         const VConjoint::EStats &es,
//                         const VConjoint::SStats &ss,
//                         const VConjoint::Options &o)
//{
//}

void VConjoint::report(const FileName &endo, const Options &o)
{
//    const auto es = analyzeE(endo, o);
//    const auto ss = analyzeS(seqs, o);
//
//    o.info("TP: " + std::to_string(ss.oc.tp()));
//    o.info("FP: " + std::to_string(ss.oc.fp()));
//    o.info("FN: " + std::to_string(ss.oc.fn()));
//
//    o.info("Generating statistics");
//
//    /*
//     * Generating VarDetect_sequins.csv
//     */
//    
//    writeQuins("VarDetect_sequins.csv", ss, o);
//
//    /*
//     * Generating VarDetect_summary.stats
//     */
//    
//    writeSummary("VarDetect_summary.stats", endo, seqs, es, ss, o);
//    
//    /*
//     * Generating VarDetect_detected.csv
//     */
//    
//    writeDetected("VarDetect_detected.csv", ss, o);
//    
//    /*
//     * Generating VarDetect_ROC.R
//     */
//    
//    o.generate("VarDetect_ROC.R");
//    o.writer->open("VarDetect_ROC.R");
//    
//    o.writer->write(createVGROC("VarDetect_detected.csv", "data$Depth", "'FP'"));
//    o.writer->close();
}
