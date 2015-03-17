#include <fstream>
#include <assert.h>
#include <algorithm>
#include "ParserFA.hpp"
#include "ParserGTF.hpp"
#include "StandardFactory.hpp"

using namespace std;

Chromosome StandardFactory::reference()
{
	std::ifstream in("/Users/user1/Sources/QA/Data/Standards/ChrT.5.10.fa");
	//std::ifstream in("C://Sources//QA//Data//Standards//ChrT.5.10.fa");
	std::string line;

	// Assume that the first line contains only the name of the chromosome
	std::getline(in, line);

    Chromosome c;
    
	// Remove the '<' prefix
	c.id = line.substr(1, line.size());
    
	c.end = std::numeric_limits<Locus>::min();
	c.start = std::numeric_limits<Locus>::max();

    ParserGTF::parse("/Users/user1/Sources/ABCD/standards/RNAstandards.gtf", [&](const Feature &f)
	//ParserGTF::parse("C://Sources//QA//Data//Standards//RNAstandards.gtf", [&](const Feature &f)
	{
		c.end = std::max(c.end, f.end);
		c.start = std::min(c.start, f.start);
        c.fs.push_back(f);
    });

    /*
     * Extract junctions. Only possible after the whole file has been parsed because the fearures might not be sorted.
     */
    
    
    
    
	ParserFA::parse("/Users/user1/Sources/QA/Data/Standards/RNAsequins.fa", [&](const Sequence &s)
	//ParserFA::parse("C://Sources//QA//Data//Standards//RNAsequins.fa", [&](const Sequence &s)
	{
		c.sequins[s.id] = s;
	});

    return c;
}