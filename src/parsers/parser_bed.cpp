#include <assert.h>
#include "data/reader.hpp"
#include "parsers/parser_bed.hpp"
#include <boost/algorithm/string.hpp>

using namespace Anaquin;

void ParserBed::parse(const Reader &r, Callback x)
{
    Data f;
    ParserProgress p;

    std::vector<std::string> sizes, starts, tokens;

    while (r.nextTokens(tokens, "\t"))
    {
        // Empty line?
        if (tokens.size() == 1)
        {
            return;
        }

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
        }

        x(f, p);
        p.i++;
    }
}