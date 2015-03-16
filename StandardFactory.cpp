#include <fstream>
#include <assert.h>
#include "ParserFA.hpp"
#include "ParserGTF.hpp"
#include "StandardFactory.hpp"

using namespace std;

Chromosome StandardFactory::reference()
{
	std::ifstream in("/Users/user1/Sources/QA/Data/Standards/ChrT.5.10.fa");
	std::string line;

	// Assume that the first line contains only the name of the chromosome
	std::getline(in, line);

    Chromosome c;
    
	// Remove the '<' prefix
	c.id = line.substr(1, line.size());
    
    ParserGTF::parse("/Users/user1/Sources/ABCD/standards/RNAstandards.gtf", [&](const Feature &f)
    {
        c.fs.push_back(f);
    });

    return c;
}

SequinMap StandardFactory::sequins()
{
    SequinMap sm;
    
    ParserFA::parse("/Users/user1/Sources/QA/Data/Standards/RNAsequins.fa", [&](const Sequence &s)
    {
        sm[s.id] = s;
    });

    return sm;
}

std::shared_ptr<Sequence> StandardFactory::sequence()
{
	std::shared_ptr<Sequence> seq;

	/*
	 * Extract the name of the in-sillico chromosome. Assume only a single chromosome in the file.
	 */

	ParserFA::parse("/Users/user1/Sources/QA/Data/Standards/ChrT.5.10.fa", [&](const Sequence &s)
	{
		seq = std::shared_ptr<Sequence>(new Sequence());
	});

	assert(seq);
	return seq;
}