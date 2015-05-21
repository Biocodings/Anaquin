#include <map>
#include "reader.hpp"
#include "tokens.hpp"
#include "parsers/parser_tmap.hpp"

using namespace Spike;

enum TMapField
{
    F_RefGeneID   = 0,
    F_RefID       = 1,
    F_ClassCode   = 2,
    F_CuffGeneID  = 3,
    F_CuffID      = 4,
    F_FMI         = 5,
    F_FPKM        = 6,
    F_FPKM_Lo     = 7,
    F_FPKM_Hi     = 8,
    F_Cov         = 9,
    F_Len         = 10,
    F_MajorIsoId  = 11,
    F_RefMatchLen = 12,
};

void ParserTMap::parse(const std::string &file, std::function<void (const TMap &, const ParserProgress &)> f)
{
    Reader i(file);

    TMap t;
    ParserProgress p;

    std::string line;
    std::vector<std::string> tokens;
    
    while (i.nextLine(line))
    {
        if (p.i++ == 0)
        {
            continue;
        }
        
        Tokens::split(line, "\t", tokens);

        t.refGenID = tokens[F_RefGeneID];
        t.refID = tokens[F_RefID];
        t.fpkm  = stof(tokens[F_FPKM]);
        t.lFPKM = stof(tokens[F_FPKM_Lo]);
        t.uFPKM = stof(tokens[F_FPKM_Hi]);

        f(t, p);
    }
}