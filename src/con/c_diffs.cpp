#include "con/c_diffs.hpp"
#include "con/c_single.hpp"
#include <ss/regression/lm.hpp>

using namespace Spike;

CDiffs::Stats CDiffs::analyze(const std::string &fileA, const std::string &fileB, const Options &options)
{
    CDiffs::Stats stats;

    /*
     * Let's reuse the code for single mixture. We'll create create a histogram for both mixtures.
     */

    options.terminal->write("Analyzing mixuture A: " + fileA);
    const auto a = CSingle::analyze(fileA);

    options.terminal->write("Analyzing mixuture B: " + fileB);
    const auto b = CSingle::analyze(fileB);

    options.terminal->write("Mergin mixtures");
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
            options.terminal->write((boost::format("Warning: %1% defined in mixture A but not in mixture B") % id).str());
            continue;
        }

        const auto baseID = s.c_map.at(id);
        
        // Calculate known fold change between mixture A and B
        const auto known = (b.expect.at(id) / a.expect.at(id));

        // Calculate actual fold change between mixture A and B
        const auto actual = (b.actual.at(id) / a.actual.at(id));

        stats.z.push_back(id);
        stats.x.push_back(log(known));
        stats.y.push_back(log(actual));
    }
    
    // Perform a linear regreession
    stats.linear();

    options.writer->open("conjoint_diffs.R");
    options.writer->write(RWriter::write(stats.x, stats.y, stats.z, "?", 0.0));
    options.writer->close();
    
	return stats;
}