#include <set>
#include <iostream>
#include "ParserGTF.hpp"
#include "AssemblyAnalyst.hpp"
#include "StandardFactory.hpp"

using namespace std;

AssemblyStats AssemblyAnalyst::analyze(const std::string &file, Sequins s, Reads n)
{
	AssemblyStats stats;
	const auto r = StandardFactory::reference();

	auto assign = [&](const Feature &f, ConfusionMatrix &m)
	{
		if (f.chromo == r.id)
		{
			if (r.loc.contains(f.loc))
			{
				if (r.matchFeature(f))
				{
					m.tp++;
				}
				else
				{
					m.fp++;
				}
			}
			else
			{
				m.fp++;
			}
		}
		else
		{
			if (r.loc.contains(f.loc))
			{
				m.fn++;
			}
			else
			{
				m.tn++;
			}
		}
	};

	ParserGTF::parse(file, [&](const Feature &f)
	{
		assign(f, stats.base);

		switch (f.type)
		{
			case Exon:
			{
				assign(f, stats.exon);
				break;
			}

			case Intron:
			{
				assign(f, stats.intron);
				break;
			}

			default: { break; }
		}
	});

	return stats;
}