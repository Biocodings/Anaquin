//#include "LadQuin/l_diffs.hpp"
//#include "LadQuin/l_express.hpp"
//
//using namespace Anaquin;
//
//LDiffs::Stats LDiffs::report(const FileName &fileA, const FileName &fileB, const Options &o)
//{
//    LDiffs::Stats stats;
//
//    // Copy the pointers across
//    auto opt = LExpress::Options();
//    
//    opt.writer = o.writer;
//    opt.logger = o.logger;
//    opt.output = o.output;
//
//    opt.mix = Mix_1;
//    const auto a = LExpress::analyze(fileA, opt);
//
//    opt.mix = Mix_2;
//    const auto b = LExpress::analyze(fileB, opt);
//
//    /*
//     * Print a warning message for each sequin detected in B but not in A
//     */
//
//    o.info("Checking for sequins in mix B but not in mix A");
//
//    for (const auto &i : b.data.normalized)
//    {
//        const auto &id = i.first;
//        
//        if (!a.data.normalized.at(id))
//        {
//            o.warn((boost::format("Warning: %1% defined in mixture B but not in mixture A") % id).str());
//        }
//    }
//    
//    o.info((boost::format("%1% sequins in mix A") % a.data.normalized.size()).str());
//    o.info((boost::format("%1% sequins in mix B") % b.data.normalized.size()).str());
//
//    /*
//     * Try for each detected sequin. But only if it's detected in both mixtures.
//     */
//    
//    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";
//
//    o.writer->open("LadderDifferent_quins.csv");
//    o.writer->write((boost::format(format) % "ID"
//                                           % "Expected A"
//                                           % "Expected B"
//                                           % "Expected B/A"
//                                           % "Measured A"
//                                           % "Measured B"
//                                           % "Measured B/A"
//                                           % "Normalized A"
//                                           % "Normalized B"
//                                           % "Normalized B/A"
//                                           % "Adjusted A"
//                                           % "Adjusted B"
//                                           % "Adjusted B/A").str());
//
//    for (const auto &i : a.data.normalized)
//    {
//        // Eg: C_02_C
//        const auto &seqID = i.first;
//
//        // Don't bother unless the sequin is detected in both mixtures
//        if (!b.data.normalized.at(seqID))
//        {
//            o.writer->write((boost::format(format) % seqID
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA"
//                                                   % "NA").str());
//            continue;
//        }
//
//        // Known fold change between mixture A and B
//        const auto known = (b.data.expect.at(seqID) / a.data.expect.at(seqID));
//
//        // Measured fold change between mixture A and B
//        const auto measured = (b.data.measured.at(seqID) / a.data.measured.at(seqID));
//
//        // Normalized fold change between mixture A and B
//        const auto normalized = (b.data.normalized.at(seqID) / a.data.normalized.at(seqID));
//
//        // Measured adjusted fold change between mixture A and B
//        const auto adjusted = (b.data.adjusted.at(seqID) / a.data.adjusted.at(seqID));
//
//        assert(a.data.measured.count(seqID) && b.data.measured.count(seqID));
//
//        stats.data[ChrT].add(seqID, known, normalized);
//        //stats.add(seqID, known, adjusted);
//
//        o.writer->write((boost::format(format) % seqID
//                                               % a.data.expect.at(seqID)
//                                               % b.data.expect.at(seqID)
//                                               % known
//                                               % a.data.measured.at(seqID)
//                                               % b.data.measured.at(seqID)
//                                               % measured
//                                               % a.data.normalized.at(seqID)
//                                               % b.data.normalized.at(seqID)
//                                               % normalized
//                                               % a.data.adjusted.at(seqID)
//                                               % b.data.adjusted.at(seqID)
//                                               % adjusted).str());
//    }
//    
//    o.writer->close();
//
//    /*
//     * Generating summary statistics
//     */
//    
//    o.info("Generating summary statistics");
//    //AnalyzeReporter::linear("LadderDifferent_summary.stats", fileA, fileB, a, b, stats, "sequins", o.writer);
//
//    /*
//     * Generating Bioconductor
//     */
//    
//    o.info("Generating Bioconductor");
//    ////AnalyzeReporter::scatter(stats, "", "LadderDifferent", "Expected concentration (attomol/ul)", "Measured coverage (q)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 reads)", o.writer, true, false);
//
//	return stats;
//}