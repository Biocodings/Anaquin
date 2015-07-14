#include "ladder/l_diffs.hpp"
#include "ladder/l_abund.hpp"
#include <ss/regression/lm.hpp>

using namespace Anaquin;

LDiffs::Stats LDiffs::analyze(const std::string &fileA, const std::string &fileB, const Options &options)
{
    LDiffs::Stats stats;

    /*
     * Let's reuse the code for single mixture. We'll create create a histogram for both mixtures.
     */

    options.output->write("Analyzing mixuture A: " + fileA);
    const auto a = LAbund::analyze(fileA);

    options.output->write("Analyzing mixuture B: " + fileB);
    const auto b = LAbund::analyze(fileB);

    options.output->write("Mergin mixtures");
    const auto &s = Standard::instance();

    /*
     * Try for each detected sequin in the experiment
     */
    
    for (const auto &i : a.actual)
    {
        const auto &id = i.first;

        // Don't bother if the sequin isn't detected in either mixture
        if (!a.actual.at(id) || !b.actual.at(id))
        {
            options.output->write((boost::format("Warning: %1% defined in mixture A but not in mixture B") % id).str());
            continue;
        }

        const auto baseID = s.l_map.at(id);

        // Calculate known fold change between mixture A and B
        const auto known = (b.expect.at(id) / a.expect.at(id));

        // Calculate actual fold change between mixture A and B
        const auto actual = (b.actual.at(id) / a.actual.at(id));

        stats.z.push_back(id);
        stats.x.push_back(log(known));
        stats.y.push_back(log(actual));
    }
    
    AnalyzeReporter::linear(stats, "ladder_diffs", "FPKM", options.writer);

	return stats;
}