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

    options.info("Merging mixtures");
    options.logInfo((boost::format("%1% sequins in mix A") % a.actual.size()).str());
    options.logInfo((boost::format("%1% sequins in mix B") % b.actual.size()).str());
    
    /*
     * Try for each detected sequin. But only if it's detected in both mixtures. Otherwise, the fold-
     * change is infinite.
     */
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";

    options.writer->open("ladder_hist.csv");
    options.writer->write((boost::format(format) % "ID"
                                                 % "abd_A"
                                                 % "abd_B"
                                                 % "exp_A"
                                                 % "exp_B"
                                                 % "exp_D"
                                                 % "mea_A"
                                                 % "mea_B"
                                                 % "adj_A"
                                                 % "adj_B"
                                                 % "adj_D").str());

    for (const auto &i : a.actual)
    {
        // Eg: C_02_C
        const auto &id = i.first;

        // Don't bother unless the sequin is detected in both mixtures
        if (!b.actual.at(id))
        {
            options.warn((boost::format("Warning: %1% defined in mixture A but not in mixture B") % id).str());
            continue;
        }

        const auto baseID = s.l_map.at(id);

        // Calculate known fold change between mixture A and B
        const auto known = (b.expect.at(id) / a.expect.at(id));

        // Calculate actual fold change between mixture A and B
        const auto measured = (b.adjusted.at(id) / a.adjusted.at(id));

        options.logInfo((boost::format("%1%\t%2%\t%3%") % id % known % measured).str());
        
        stats.z.push_back(id);
        stats.x.push_back(log2(known));
        stats.y.push_back(log2(measured));
        
        assert(a.measure.count(id) && b.measure.count(id));
        
        options.writer->write((boost::format(format) % id
                                                     % a.measure.at(id)
                                                     % b.measure.at(id)
                                                     % a.expect.at(id)
                                                     % b.expect.at(id)
                                                     % known
                                                     % a.actual.at(id)
                                                     % b.actual.at(id)
                                                     % a.adjusted.at(id)
                                                     % b.adjusted.at(id)
                                                     % measured).str());
    }
    
    options.writer->close();
    options.logInfo("Checking for sequins in mix B but not in mix A");

    /*
     * Print a warning message for each sequin detected in B but not in A
     */
    
    for (const auto &i : b.actual)
    {
        const auto &id = i.first;

        if (!a.actual.at(id))
        {
            const auto msg = (boost::format("Warning: %1% defined in mixture B but not in mixture A") % id).str();
            options.out(msg);
            options.warn(msg);
        }
    }

    options.info("Generating linear model");
    AnalyzeReporter::linear(stats, "ladder_diffs", "FPKM", options.writer);
    
	return stats;
}