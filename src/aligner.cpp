#include <map>
#include <assert.h> 
#include "types.hpp"
#include "aligner.hpp"
#include "reader_fa.hpp"
#include "sam_reader.hpp"

using namespace std;

AlignerStatistics Aligner::analyze(const std::string &file)
{
	AlignerStatistics stats;

	Reads tp = 0;
	Reads tn = 0;
	Reads fp = 0;
	Reads fn = 0;

	std::map<ID, Sequence> seqs;

	ReaderFA::read("C:\\Sources\\QA\\data\\sequins\\RNAsequins.fa", [&](const Sequence &s)
	{
		seqs[s.id] = s;
	});

	/*
	 * For each alignment, we'd want to know if it's aligned with a sequin. If it's, check whether this is really
	 * a correct alignment. If it's not aligned with a sequin, we'd still need to check if it's supposed to be
	 * aligned.
	 */

	SAMReader::read(file, [&](const Alignment &align)
	{
		if (align.mapped)
		{
			// Aligned to a sequin
			if (align.name == "chrT")
			{
				/*
				 * Is this a correct alignment? Assume it's for now... Most likely need to compare with the raw DNA sequence!
				 */

				tp++;

				// if not -> fp++

				stats.n_si++;
			}

			// Aligned to a real sample
			else
			{
				// if really sequin then fn++
				// if really NOT sequin then tn++

				stats.n_sa++;
			}

			stats.n++;
		}
	});

	stats.p_si = static_cast<Percentage>(stats.n_si / stats.n);
	stats.p_sa = static_cast<Percentage>(stats.n_sa / stats.n);
	stats.dilution = static_cast<Percentage>(stats.n_si / stats.n_sa);

	return stats;
}
