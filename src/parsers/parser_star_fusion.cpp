#include <stdexcept>
#include "data/tokens.hpp"
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

enum Fields
{
    FusionName,
    JunctionReads,
    SpanningFrags,
    LeftGene,
    LeftBreakpoint,
    LeftDistFromRefExonSplice,
    RightGene,
    RightBreakpoint,
    RightDistFromRefExonSplice,
};

static void parseBreak(const std::string &s, std::string &chr, Base &l, Strand &o)
{
    /*
     * Eg: chrT:9035684:-
     */
    
    std::vector<std::string> tokens;
    Tokens::split(s, ":", tokens);

    assert(tokens.size() == 3);
    
    // Eg: chrT
    chr = tokens[0];

    // Eg: 9035684
    l = std::stoi(tokens[1]);
    
    // Eg: '-'
    o = tokens[2] == "-" ? Backward : Forward;
}

void ParserStarFusion::parse(const Reader &r, Functor f)
{
    std::string line;
    std::vector<std::string> tokens;

    Data data;
    ParserProgress p;
    
    while (r.nextLine(line))
    {
        if (line[0] == '#')
        {
            continue;
        }

        p.i++;
        Tokens::split(line, "\t", tokens);
        
        const auto leftGene  = tokens[LeftGene];
        const auto rightGene = tokens[RightGene];

        parseBreak(tokens[LeftBreakpoint],  data.cID_1, data.l1, data.s1);
        parseBreak(tokens[RightBreakpoint], data.cID_2, data.l2, data.s2);

        // Measured abundance
        data.reads = stoi(tokens[JunctionReads]);

        f(data, p);
    }
}