#include "ParserGTF.hpp"
#include "SillicoFactory.hpp"
#include "AssemblyAnalyst.hpp"

using namespace std;

AssemblyStats AssemblyAnalyst::analyze(const std::string &file)
{
	AssemblyStats stats;
	//const auto sillico = SillicoFactory::sequence();

	/*
	 * Read for the assembled transcript
	 */

	struct TranscriptReader : public ReaderGTF
	{
		void exon(const Feature &f) override
		{

		}

		AssemblyStats *stats;
	};

	TranscriptReader tReader;
	ParserGTF::parse(SillicoFactory::transGTF(), tReader);

	/*
	 * Check for the features in the in-sillico chromosome, have they been assembled?
	 */

	struct SillicoReader : public ReaderGTF
	{
		void exon(const Feature &f) override
		{
			/*
			 * Check if the known exon has been assembled
			 */
		}

		AssemblyStats *stats;
	};

	SillicoReader sReader;
	sReader.stats = &stats;

	// Check the features listed in the sillico transcripts
	ParserGTF::parse(SillicoFactory::transGTF(), sReader);

	return stats;
}