#include "tools/tools.hpp"
#include "VarQuin/v_conjoin.hpp"

using namespace Anaquin;

extern Scripts PlotVLODR();
extern Scripts PlotVGROC();
extern Scripts PlotVCROC();
extern Scripts PlotAllele();

extern Path __output__;
extern std::string __full_command__;

VConjoint::Stats VConjoint::analyze(const FileName &file, const Options &o)
{
//    const auto &r = Standard::instance().r_var;

    VConjoint::Stats stats;
    
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

void VConjoint::report(const FileName &file, const Options &o)
{
//    const auto stats = analyze(file, o);
//
//    o.info("Generating statistics");
//
//    /*
//     * Generating VarConjoin_sequins.csv
//     */
//    
//    writeQuins("VarConjoin_sequins.csv", ss, o);
//
//    /*
//     * Generating VarConjoin_summary.stats
//     */
//    
//    writeSummary("VarConjoin_summary.stats", endo, seqs, es, ss, o);
//    
//    /*
//     * Generating VarConjoin_detected.csv
//     */
//    
//    writeDetected("VarConjoin_detected.csv", ss, o);
//    
//    /*
//     * Generating VarConjoin_ROC.R
//     */
//    
//    o.generate("VarConjoin_ROC.R");
//    o.writer->open("VarConjoin_ROC.R");
//    
//    o.writer->write(createVGROC("VarConjoin_detected.csv", "data$Depth", "'FP'"));
//    o.writer->close();
}
