#include <assert.h>
#include "ParserFA.hpp"
#include "SillicoFactory.hpp"

using namespace std;

std::string SillicoFactory::transGTF()
{
	return "C:\\Sources\\QA\\data\\standards\\RNAstandards.gtf";
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