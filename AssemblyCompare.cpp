#include "ParserGTF.hpp"
#include "SillicoFactory.hpp"
#include "AssemblyCompare.hpp"

using namespace std;

AssemblyStatistics AssemblyCompare::analyze(const std::string &file)
{
	AssemblyStatistics stats;
	const auto sillico = SillicoFactory::sequence();

	/*
	 * Construct a data-structure for the assembled transcript
	 */

	struct AssemblyReader : public ReaderGTF
	{
		void exon(const Feature &f) override
		{

		}

		AssemblyStatistics *stats;
	};

	/*
	 * Check for the features in the in-sillico chromosome, have they been assembled?
	 */

	struct AssemblyReader : public ReaderGTF
	{
		void exon(const Feature &f) override
		{
			/*
			 * Check if the known exon has been assembled
			 */
		}

		AssemblyStatistics *stats;
	};

	AssemblyReader reader;
	reader.stats = &stats;

	// Check the features listed in the sillico transcripts
	ParserGTF::parse(SillicoFactory::transGTF(), reader);

	return stats;
}