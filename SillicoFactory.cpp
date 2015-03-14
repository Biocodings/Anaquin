#include <assert.h>
#include "ParserFA.hpp"
#include "ParserGTF.hpp"
#include "SillicoFactory.hpp"

using namespace std;

std::string SillicoFactory::transGTF()
{
	return "/Users/user1/Sources/ABCD/standards/RNAstandards.gtf";
}

std::shared_ptr<Sequence> SillicoFactory::sequence()
{
	std::shared_ptr<Sequence> seq;

	/*
	* Extract the name of the in-sillico chromosome. Assume only a single chromosome in the file.
	*/

	ParserFA::parse("C:\\Sources\\QA\\data\\standards\\ChrT.5.10.fa", [&](const Sequence &s)
	{
		seq = std::shared_ptr<Sequence>(new Sequence());
	});

	assert(seq);
	return seq;
}

FeatureMap SillicoFactory::features()
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





