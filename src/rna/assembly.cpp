#include <limits>
#include <iostream>
#include "assembly.hpp"
#include "statistics.hpp"
#include "standard_factory.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Spike;

AssemblyStats Assembly::analyze(const std::string &file)
{
    AssemblyStats stats;
    const auto r = StandardFactory::reference();

	ParserGTF::parse(file, [&](const Feature &f, ParserProgress &p)
	{
        if (false /*p.i > n*/)
        {
        }
        else
        {
            // Binary classification for the base-level
            classify(r.fs, r, f, stats.base);
            
            switch (f.type)
            {
                case Exon: { classify(r.fs, r, f, stats.exon);   break; }
                default:   { break; }
            }
        }
	});

	return stats;
}
