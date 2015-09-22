#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "parsers/parser_blast.hpp"

using namespace Anaquin;

enum PSLField
{
    PSL_Matches,
    PSL_MisMatches,
    PSL_RepMatches,
    PSL_NCount,
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

void ParserBlast::parse(const Reader &r, Callback c)
{
    ParserProgress p;

    std::string line;
    std::vector<std::string> toks;

    BlastLine l;
    
    while (r.nextLine(line))
    {
        if (p.i++ <= 3)
        {
            continue;
        }
        
        Tokens::split(line, "\t", toks);

        l.qName     = toks[PSL_QName];
        l.tName     = toks[PSL_TName];

        l.tStart    = stoi(toks[PSL_TStart]);
        l.tEnd      = stoi(toks[PSL_TEnd]);
        l.tSize     = stoi(toks[PSL_TSize]);

        l.qStart    = stoi(toks[PSL_QStart]);
        l.qEnd      = stoi(toks[PSL_QEnd]);
        l.qSize     = stoi(toks[PSL_QSize]);

        l.match     = stoi(toks[PSL_Matches]);
        l.mismatch  = stoi(toks[PSL_MisMatches]);
        
        l.qGap      = stoi(toks[PSL_QGap_Bases]);
        l.tGap      = stoi(toks[PSL_TGap_Bases]);

        l.qGapCount = stoi(toks[PSL_QGap_Count]);
        l.tGapCount = stoi(toks[PSL_TGap_Count]);

        c(l, p);
    }
}