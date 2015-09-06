#include "variant/v_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

VAlign::Stats VAlign::analyze(const std::string &file, const Options &options)
{
    VAlign::Stats stats;
//    static const auto &s = Standard::instance();
//
//    options.info("Parsing alignment file");
//
//    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &p)
//    {
//        if (!align.i && (p.i % 1000000) == 0)
//        {
//            options.wait(std::to_string(p.i));
//        }
//        
//        if (align.id != s.id && !align.i)
//        {
//            stats.n_genome++;
//        }
//        
//        if (!align.mapped || align.id != s.id)
//        {
//            return;
//        }
//        else if (!align.i)
//        {
//            stats.n_chrT++;
//        }
//
//        /*
//         * Collect statistics for the alignment
//         */
//        
//        Feature f;
//
//        if (classify(stats.p.m, align, [&](const Alignment &)
//        {
//            return find(s.fs_1.begin(), s.fs_1.end(), align, f) ? Positive : Negative;
//        }))
//        {
//            stats.h_seq.at(f.tID)++;
//            stats.h_base.at(s.seq2base.at(f.tID))++;
//        }
//    });
//
//    sums(stats.h_seq, stats.p.m.nr);
//
//    /*
//     * Generate an abundance plot for the accuracy of quantification, the measured DNA
//     * standard abundance (in FPKM) relative to the known concentration (in attamoles/ul)
//     * of each DNA standard.
//     */
//
//    // The total number of reads aligned
////    const auto n = stats.n_chrT;
//
//    // Known concentration for the given mixture
////    const auto &m = s.v_seq(options.mix);
////
////    for (const auto &i : stats.c)
////    {
////        if (!i.second)
////        {
////            continue;
////        }
////        
////        const auto s = m.at(i.first);
////        
////        // Compare the FPKM with the known concentration
////        const auto known = s.abund();
////        
////        // Calculate FPKM for the sequin
////        const double measured = (std::pow(10, 9) * static_cast<double>(i.second)) / (n * s.l.length());
////
////        stats.z.push_back(i.first);
////        stats.x.push_back(log2(known));
////        stats.y.push_back(log2(measured));
////    }
////    
////    // Perform a linear regreession
////    stats.linear();
////
////    options.info("Generating statistics");
////    
//    // Calculate for the sensitivity
//    stats.p.s = Expression::analyze(stats.h_seq, s.seqs_1);
//
//    //AnalyzeReporter::stats("VarAlign_summary.stats", stats.p, stats.h_seq, options.writer);
//
	return stats;
}
