#include "ladder/l_diffs.hpp"
#include "ladder/l_abund.hpp"
#include <ss/regression/lm.hpp>

using namespace Anaquin;

LDiffs::Stats LDiffs::analyze(const std::string &fileA, const std::string &fileB, const Options &options)
{
    LDiffs::Stats stats;

    // Copy the pointers across
    auto opt = LAbund::Options();
    
    opt.writer = options.writer;
    opt.logger = options.logger;
    opt.output = options.output;

    options.info("Analyzing mixuture A: " + fileA);
    const auto a = LAbund::analyze(fileA, opt);

    opt.mix = MixB;
    options.info("Analyzing mixuture B: " + fileB);
    const auto b = LAbund::analyze(fileB, opt);

    const auto &s = Standard::instance();

    options.logInfo("Checking for sequins in mix B but not in mix A");
    
    /*
     * Print a warning message for each sequin detected in B but not in A
     */
    
    for (const auto &i : b.normalized)
    {
        const auto &id = i.first;
        
        if (!a.normalized.at(id))
        {
            const auto msg = (boost::format("Warning: %1% defined in mixture B but not in mixture A") % id).str();
            options.out(msg);
            options.warn(msg);
        }
    }
    
    options.info("Merging mixtures");
    options.logInfo((boost::format("%1% sequins in mix A") % a.normalized.size()).str());
    options.logInfo((boost::format("%1% sequins in mix B") % b.normalized.size()).str());
    
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
        const auto &id = i.first;

        // Don't bother unless the sequin is detected in both mixtures
        if (!b.normalized.at(id))
        {
            options.warn((boost::format("Warning: %1% defined in mixture A but not in mixture B") % id).str());
            continue;
        }

        const auto baseID = s.l_map.at(id);

        // Calculate known fold change between mixture A and B
        const auto known = (b.expect.at(id) / a.expect.at(id));

        // Calculate actual fold change between mixture A and B
        const auto adjusted = (b.adjusted.at(id) / a.adjusted.at(id));

        options.logInfo((boost::format("%1%\t%2%\t%3%") % id % known % adjusted).str());

        stats.z.push_back(id);
        stats.x.push_back(log2(known));
        stats.y.push_back(log2(adjusted));
        
        assert(a.measured.count(id) && b.measured.count(id));
        
        options.writer->write((boost::format(format) % id
                                                     % a.expect.at(id)
                                                     % b.expect.at(id)
                                                     % known
                                                     % a.measured.at(id)
                                                     % b.measured.at(id)
                                                     % a.normalized.at(id)
                                                     % b.normalized.at(id)
                                                     % a.adjusted.at(id)
                                                     % b.adjusted.at(id)
                                                     % adjusted).str());
    }
    
    options.writer->close();

    AnalyzeReporter::linear(stats, "LadderDifferent", "FPKM", options.writer, true, true, false);

	return stats;
}