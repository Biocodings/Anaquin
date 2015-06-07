#include <assert.h>
#include "data/reader.hpp"
#include "parsers/parser_bed.hpp"
#include <boost/algorithm/string.hpp>

using namespace Spike;

void ParserBED::parse(const Reader &r, Callback x)
{
    BedFeature bf;
    ParserProgress p;

    /*
     *     Required
     *    ----------
     *
     *    1. Name of the chromosome
     *    2. Starting position
     *    3. Ending position
     *
     *     Optional
     *    ----------
     *
     *    4. Name of the BED line
     *    5. Score
     *    6. Strand
     *    7. Thick starting position
     *    8. Thick ending position
     *    9. An RGB value
     *   10. The number of blocks (exons) in the BED line
     *   11. A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
     *   12. A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart.
     *       The number of items in this list should correspond to blockCount.
     */

    std::vector<std::string> sizes, starts, tokens;

    while (r.nextTokens(tokens, "\t"))
    {
        // Empty line?
        if (tokens.size() == 1)
        {
            return;
        }
        
        bf.blocks.clear();
        
        // Name of the chromosome
        bf.id = tokens[0];

        // Position of the feature in standard chromosomal coordinates
        bf.l  = Locus(stod(tokens[1]) + 1, stod(tokens[2]));

        if (tokens.size() >= 4)
        {
            // Name of the BED line (eg: gene)
            bf.name = tokens[3];

            if (tokens.size() >= 10)
            {
                boost::split(sizes,  tokens[10], boost::is_any_of(","));
                boost::split(starts, tokens[11], boost::is_any_of(","));
                assert(sizes.size() == starts.size());

                for (auto i = 0; i < stod(tokens[9]); i++)
                {
                    const BasePair start = stod(starts[i]);
                    const BasePair size  = stod(sizes[i]);

                    bf.blocks.push_back(Locus(bf.l.start + start, bf.l.start + start + size - 1));
                }
            }
        }

        x(bf, p);
        p.i++;
    }
}