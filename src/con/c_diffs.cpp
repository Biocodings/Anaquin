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
        if (!b.s_correct.count(id))
        {
            std::cout << "Warning: " + i.first << std::endl;
            continue;
        }
        
        const auto &id = i.first;
        
        assert(b.s_correct.count(id));

        // Calculate actual fold change between mixture A and B
        const auto actual = log(b.s_correct.at(id) / a.s_correct.at(id));
        
        // Calculate known fold change between mixture A and B
        const auto known = log(s.c_seqs_B.at(id).abund() / s.c_seqs_A.at(id).abund());

        stats.x.push_back(known);
        stats.y.push_back(actual);
        stats.z.push_back(id);
    }
    
    // Perform a linear regreession
    stats.linear();

    options.writer->open("conjoint_diffs.R");
    options.writer->write(RWriter::write(stats.x, stats.y, stats.z, "?", 0.0));
    options.writer->close();
    
	return stats;
}