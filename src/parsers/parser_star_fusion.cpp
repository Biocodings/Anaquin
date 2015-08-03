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

static void parseBreak(const std::string &s, std::string &chr, BasePair &l, Strand &o)
{
    /*
     * Eg: chrT:9035684:-
     */
    
    std::vector<std::string> tokens;
    Tokens::split(s, ":", tokens);

    assert(tokens.size() == 3);
    
    // Eg: chrT
    chr = tokens[0];
    
    l = std::stoi(tokens[1]);
    
    // Eg: '-'
    o = tokens[2] == "-" ? Backward : Forward;
}

void ParserStarFusion::parse(const Reader &r, Functor f)
{
    std::string line;
    std::vector<std::string> tokens;

    Fusion d;
    ParserProgress p;
    
    while (r.nextLine(line))
    {
        if (line[0] == '#')
        {
            continue;
        }
        
        p.i++;
        Tokens::split(line, "\t", tokens);
        
        const auto leftGene  = tokens[3];
        const auto rightGene = tokens[6];

        parseBreak(tokens[4], d.l_chr, d.l_break, d.l_strand);
        parseBreak(tokens[7], d.r_chr, d.r_break, d.r_strand);

        f(d, p);
    }
}