#include "conj/c_join.hpp"
#include "stats/expression.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

CJoin::Stats CJoin::analyze(const std::string &file, const Options &options)
{
    CJoin::Stats stats;

    // Construct a histogram of the aligned sequin
    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        if (align.i == 0)
        {
            stats.hist[align.id]++;
        }
    });

    options.writer->open("conjoin_histogram.stats");
    
    for (const auto &i : stats.hist)
    {
        options.writer->write((boost::format("%1%\t%2%") % i.first % i.second).str());
    }

    options.writer->close();
    
	return stats;
}