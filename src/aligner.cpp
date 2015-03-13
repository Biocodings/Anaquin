#include "types.hpp"
#include "aligner.hpp"
#include "sam_reader.hpp"

using namespace std;

AlignerStatistics Aligner::analyze(const std::string &file)
{
	AlignerStatistics stats;

	SAMReader::read(file, [&](const Alignment &align)
	{
		if (align.name == "chrT")
		{
			stats.n_si++;
		}
		else
		{
			stats.n_sa++;
		}

		stats.n++;
	});

	stats.p_si = stats.n_si / stats.n;
	stats.p_sa = stats.n_sa / stats.n;
	stats.dilution = static_cast<Percentage>(stats.n_si / stats.n_sa);

	return stats;
}