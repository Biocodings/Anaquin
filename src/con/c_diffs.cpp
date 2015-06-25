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

    options.terminal->write("Analyzing mixuture A: " + fileA);
    const auto a = CSingle::analyze(fileA);

    options.terminal->write("Analyzing mixuture A: " + fileB);
    const auto b = CSingle::analyze(fileB);
    
    const auto &s = Standard::instance();
    int j = 0;
    for (const auto &i : a.s_correct)
    {
        std::cout << i.first << std::endl;
        const auto &id = i.first;
        
        assert(b.s_correct.count(id));

        // Calculate actual fold change between mixture A and B
        const auto actual = log((j++ + 10) + b.s_correct.at(id) / a.s_correct.at(id));
        
        // Calculate known fold change between mixture A and B
        const auto known = log((j++ + 10) + s.c_seqs_B.at(id).abund() / s.c_seqs_A.at(id).abund());

        stats.x.push_back(j++);
        stats.y.push_back(j++);
        stats.z.push_back(id);
    }

    stats.x.push_back(j++);
    stats.y.push_back(j++);
    stats.x.push_back(j++);
    stats.y.push_back(j++);
    
    // Perform a linear regreession
    stats.linear();

    options.writer->open("conjoint_diffs.R");
    options.writer->write(RWriter::write(stats.x, stats.y, stats.z, "?", 0.0));
    options.writer->close();
    
	return stats;
}