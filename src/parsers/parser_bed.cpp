#include <assert.h>
#include "data/reader.hpp"
#include "parsers/parser_bed.hpp"
#include <boost/algorithm/string.hpp>

using namespace Anaquin;

void ParserBED::parse(const Reader &r, Callback x)
{
    ParserBED::Annotation f;
    ParserProgress p;

    std::vector<std::string> sizes, starts, tokens;

    while (r.nextTokens(tokens, "\t"))
    {
        // Empty line?
        if (tokens.size() == 1)
        {
            return;
        }

        f.blocks.clear();

        // Name of the chromosome
        f.id = tokens[0];

        // Position of the feature in standard chromosomal coordinates
        f.l = Locus(stod(tokens[1]) + 1, stod(tokens[2]));

        if (tokens.size() >= 6)
        {
            // Defines the strand, either '+' or '-'
            f.strand = tokens[5] == "+" ? Forward : Backward;
        }
        
        if (tokens.size() >= 4)
        {
            // Name of the BED line (eg: gene)
            f.name = tokens[3];

            if (tokens.size() >= 10)
            {
                boost::split(sizes,  tokens[10], boost::is_any_of(","));
                boost::split(starts, tokens[11], boost::is_any_of(","));
                assert(sizes.size() == starts.size());

                for (auto i = 0; i < stod(tokens[9]); i++)
                {
                    const Base start = stod(starts[i]);
                    const Base size  = stod(sizes[i]);

                    f.blocks.push_back(Locus(f.l.start + start, f.l.start + start + size - 1));
                }
            }
        }

        x(f, p);
        p.i++;
    }
}