#include "variant/variant.hpp"
#include "variant/v_allele.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

VAllele::Stats VAllele::analyze(const FileName &file, const Options &o)
{
    VAllele::Stats stats;
    const auto &r = Standard::instance().r_var;
    
    ParserVCF::parse(file, [&](const VCFVariant &v, const ParserProgress &)
    {
        if (v.id == ChrT)
        {
            // Known coverage for allele frequnece
            const auto known = r.alleleFreq(Mix_1, match->bID);

            // Measured coverage is the number of base calls aligned and used in variant calling
            const auto measured = static_cast<double>(v.dp_a) / (v.dp_r + v.dp_a);
            
            
            
        }
    });
    
    
//    
//    parseVCF(file, stats, [&](const VCFVariant &v, const Variation *match)
//    {
//
//        /*
//         * Plotting the relative allele frequency that is established by differences
//         * in the concentration of reference and variant DNA standards.
//         */
//        
//        // Eg: D_1_12_R_373892_G/A
//        const auto id = (boost::format("%1%_%2%_%3%_%4%:") % match->id
//                                                           % match->ref
//                                                           % match->l.start
//                                                           % match->alt).str();
//        
//        //stats.chrT->h.at(match->id)++;
//        
//        // TODO: How to handle a case where variant is reported but with zero counts?
//        if (v.dp_a == 0)
//        {
//            return;
//        }
//
//        stats.data.add(id, known, measured);
//    });
    
    //stats.data.ss = r.limit(stats.chrT->h);
    stats.data.sn = static_cast<double>(stats.data.detected) / r.countVars();
    
    return stats;
}

void VAllele::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file, o);
    
//    o.writer->open("VarAllele_false.stats");
//
//    classify(file, stats, [&](const VCFVariant &v, const Variation *match)
//    {
//        // The known coverage for allele frequnece
//        const auto known = r.alleleFreq(Mix_1, match->bID);
//
//        // The measured coverage is the number of base calls aligned and used in variant calling
//        const auto measured = static_cast<double>(v.dp_a) / (v.dp_r + v.dp_a);
//
//        /*
//         * Plotting the relative allele frequency that is established by differences
//         * in the concentration of reference and variant DNA standards.
//         */
//
//        // Eg: D_1_12_R_373892_G/A
//        const auto id = (boost::format("%1%_%2%_%3%_%4%:") % match->id
//                                                           % match->ref
//                                                           % match->l.start
//                                                           % match->alt).str();
//
//        stats.chrT->h.at(match->id)++;
//
//        // TODO: How to handle a case where variant is reported but with zero counts?
//        if (v.dp_a == 0)
//        {
//            return;
//        }
//
//        stats.chrT->add(id, known, measured);
//    });
// 
//    stats.chrT->ss = r.limit(stats.chrT->h);
//    stats.chrT->sn = static_cast<double>(stats.chrT->detected) / r.countVars();
//    
    /*
     * Generate summary statistics
     */

    o.info("Generating summary statistics");
    //AnalyzeReporter::linear("VarAllele_summary.stats", file, stats, "variants", o.writer);

    //AnalyzeReporter::scatter(stats, "", "VarAllele", "Expected allele frequency (proportion)", "Measured allele frequency (proportion)", "Expected allele frequency (proportion)", "Measured allele frequency (proportion)", o.writer, false);
}