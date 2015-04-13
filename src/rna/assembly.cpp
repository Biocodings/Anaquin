#include "assembly.hpp"
#include "statistics.hpp"
#include "standard_factory.hpp"
#include "parsers/parser_gtf.hpp"

using namespace Spike;

AssemblyStats Assembly::analyze(const std::string &file, const AssemblyOptions &options)
{
    AssemblyStats stats;
    const auto r = StandardFactory::reference();

	ParserGTF::parse(file, [&](const Feature &f, ParserProgress &p)
	{
        classify(r, f, stats.base, [&](const Feature &)
        {
            return find(r.fs, f);
        });

        switch (f.type)
        {
        //    case Exon: { classify(r.fs, r, f, stats.exon);   break; }
            default:   { break; }
        }
	});

	return stats;
}
