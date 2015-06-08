#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "parsers/parser_blat.hpp"

using namespace Spike;

enum PSLField
{
    PSL_Match,
    PSL_MisMatch,
    PSL_RepMatch,
    PSL_Ns,
    PSL_QGap_Count,
    PSL_QGap_Bases,
    PSL_TGap_Count,
    PSL_TGap_Bases,
    PSL_Strand,
    PSL_QName,
    PSL_QSize,
    PSL_QStart,
    PSL_QEnd,
    PSL_TName,
    PSL_TSize,
    PSL_TStart,
    PSL_TEnd,
    PSL_Block_Count,
    PSL_Block_Sizes,
    PSL_Q_Starts,
    PSL_T_Starts
};

void ParserBlat::parse(const Reader &r, Callback c)
{
    ParserProgress p;

    std::string line;
    std::vector<std::string> tokens;

    BlatLine l;
    
    while (r.nextLine(line))
    {
        if (p.i++ <= 3)
        {
            continue;
        }
        
        Tokens::split(line, "\t", tokens);

        l.qName  = tokens[PSL_QName];
        l.tName  = tokens[PSL_TName];
        l.tStart = stoi(tokens[PSL_TStart]);
        l.tEnd   = stoi(tokens[PSL_TEnd]);

        c(l, p);
    }
}


