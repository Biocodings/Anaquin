#include "con/c_diffs.hpp"
#include "con/c_single.hpp"
#include <ss/regression/lm.hpp>
#include "parsers/parser_sam.hpp"

using namespace Spike;

CDiffs::Stats CDiffs::analyze(const std::string &fileA, const std::string &fileB, const Options &options)
{
    CDiffs::Stats stats;

    /*
     * Let's reuse the code for single mixture. We'll create create a histogram for both mixtures.
     */

    options.terminal->write("Analyzing mixuture B: " + fileB);
    const auto b = CSingle::analyze(fileB);

    options.terminal->write("Analyzing mixuture A: " + fileA);
    const auto a = CSingle::analyze(fileA);

    const auto &s = Standard::instance();

    for (const auto &i : a.s_correct)
    {
        const auto &id = i.first;
        
        if (!b.s_correct.count(i.first))
        {
            options.terminal->write((boost::format("Warning: %1% defined in mixture A but not in mixture B") % id).str());
            continue;
        }

        assert(b.s_correct.count(id));

        const auto len_a = s.c_seqs_A.at(id).length;
        const auto len_b = s.c_seqs_B.at(id).length;

        // Calculate actual fold change between mixture A and B
        const auto actual = (b.abund.at(id) / len_b) / (a.abund.at(id) / len_a);
        
        // Calculate known fold change between mixture A and B
        const auto known = (s.c_seqs_B.at(id).abund() / len_b) / (s.c_seqs_A.at(id).abund() / len_a);

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