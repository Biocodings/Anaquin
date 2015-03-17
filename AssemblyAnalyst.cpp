#include <iostream>
#include "ParserGTF.hpp"
#include "Statistics.hpp"
#include "AssemblyAnalyst.hpp"
#include "StandardFactory.hpp"

AssemblyStats AssemblyAnalyst::analyze(const std::string &file, Sequins s, Reads n)
{
    AssemblyStats stats;
    
    // Reference standards
    const auto r = StandardFactory::reference();

	ParserGTF::parse(file, [&](const Feature &f, ParserProgress &p)
	{
        if (p.i > n)
        {
            p.terminate = true;
        }
        else
        {
            classify(r.fs, r, f, stats.base);
            
            switch (f.type)
            {
                case Exon:   { classify(r.fs, r, f, stats.exon);   break; }
                case Intron: { classify(r.fs, r, f, stats.intron); break; }
                default:     { break; }
            }
        }
	});

	return stats;
}