#include "ladder/l_diffs.hpp"
#include "ladder/l_abund.hpp"
#include <ss/regression/lm.hpp>

using namespace Anaquin;

LDiffs::Stats LDiffs::report(const FileName &fileA, const FileName &fileB, const Options &options)
{
    LDiffs::Stats stats;

    // Copy the pointers across
    auto opt = LAbund::Options();
    
    opt.writer = options.writer;
    opt.logger = options.logger;
    opt.output = options.output;

    options.info("Analyzing mixuture A: " + fileA);
    const auto a = LAbund::report(fileA, opt);

    //opt.mix = Mix_2; TODO: FIX THIS WILL NOT WORK
    options.info("Analyzing mixuture B: " + fileB);
    const auto b = LAbund::report(fileB, opt);

    options.logInfo("Checking for sequins in mix B but not in mix A");
    
    /*
     * Print a warning message for each sequin detected in B but not in A
     */
    
    for (const auto &i : b.normalized)
    {
        const auto &id = i.first;
        
        if (!a.normalized.at(id))
        {
            options.warn((boost::format("Warning: %1% defined in mixture B but not in mixture A") % id).str());
        }
    }
    
    options.info("Merging mixtures");
    options.info((boost::format("%1% sequins in mix A") % a.normalized.size()).str());
    options.info((boost::format("%1% sequins in mix B") % b.normalized.size()).str());

    /*
     * Try for each detected sequin. But only if it's detected in both mixtures. Otherwise, the fold-
     * change is infinite.
     */
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";

    options.writer->open("LadderDifferent_hist.csv");
    options.writer->write((boost::format(format) % "id"
                                                 % "expect_A"
                                                 % "expect_B"
                                                 % "expect_D"
                                                 % "measure_A"
                                                 % "measure_B"
                                                 % "norm_A"
                                                 % "norm_B"
                                                 % "adjust_A"
                                                 % "adjust_B"
                                                 % "adjust_D").str());

    for (const auto &i : a.normalized)
    {
        // Eg: C_02_C
        const auto &seqID = i.first;

        // Don't bother unless the sequin is detected in both mixtures
        if (!b.normalized.at(seqID))
        {
            options.warn((boost::format("Warning: %1% defined in the first sample but not in the second sample") % seqID).str());
            continue;
        }

        // Calculate known fold change between mixture A and B
        const auto known = (b.expect.at(seqID) / a.expect.at(seqID));

        // Calculate actual fold change between mixture A and B
        const auto adjusted = (b.adjusted.at(seqID) / a.adjusted.at(seqID));

        options.logInfo((boost::format("%1%\t%2%\t%3%") % seqID % known % adjusted).str());

        stats.add(seqID, log2(known), log2(adjusted));
        assert(a.measured.count(seqID) && b.measured.count(seqID));
        
        options.writer->write((boost::format(format) % seqID
                                                     % a.expect.at(seqID)
                                                     % b.expect.at(seqID)
                                                     % known
                                                     % a.measured.at(seqID)
                                                     % b.measured.at(seqID)
                                                     % a.normalized.at(seqID)
                                                     % b.normalized.at(seqID)
                                                     % a.adjusted.at(seqID)
                                                     % b.adjusted.at(seqID)
                                                     % adjusted).str());
    }
    
    options.writer->close();

    //AnalyzeReporter::linear(stats, "LadderDifferent", "FPKM", options.writer, true, true, false);

	return stats;
}