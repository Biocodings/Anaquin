#include <fstream>
#include <assert.h>
#include "ParserFA.hpp"
#include "ParserGTF.hpp"
#include "StandardFactory.hpp"

using namespace std;

std::string StandardFactory::chromoName()
{
	std::ifstream in("/Users/user1/Sources/QA/Data/Standards/ChrT.5.10.fa");
	std::string line;

	// Assume the first line contains only the name of the chromosome
	std::getline(in, line);

	// Remove the '<' prefix
	return line.substr(1, line.size());
}

std::string StandardFactory::transGTF()
{
	return "/Users/user1/Sources/ABCD/standards/RNAstandards.gtf";
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

FeatureMap StandardFactory::features()
{
    FeatureMap mapper;
    
    struct PrivateReader : public FeatureReader
    {
        void all(const Feature &f)
        {
            (*mapper)[f.pos] = f;
            (*mapper)[f.pos+1] = f; // Hack?
            (*mapper)[f.pos-1] = f; // Hack?
        }
        
        FeatureMap *mapper;
    };
    
    PrivateReader reader;
    reader.mapper = &mapper;
    ParserGTF::parse("/Users/user1/Sources/ABCD/standards/RNAstandards.gtf", reader);

    return mapper;
}





