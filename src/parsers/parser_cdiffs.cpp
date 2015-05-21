#include "reader.hpp"
#include "tokens.hpp"
#include "parsers/parser_cdiffs.hpp"
#include <boost/algorithm/string/predicate.hpp>

using namespace SS;
using namespace Spike;

enum TrackingField
{
    FTestID     = 0,
    FGeneID     = 1,
    FStatus     = 6,
    FFPKM_1     = 7,
    FFPKM_2     = 8,
    FLogFold    = 9,
    FTestStats  = 10,
    FPValue     = 11,
    FQValue     = 12,
};

void ParserCDiffs::parse(const std::string &file, std::function<void (const TrackingDiffs &, const ParserProgress &)> f)
{
    Reader i(file);
    TrackingDiffs t;
    ParserProgress p;

    std::string line;
    std::vector<std::string> tokens;
    
    while (i.nextLine(line))
    {
        p.i++;

        if (boost::starts_with(line, "test_id"))
        {
            continue;
        }
        
        Tokens::split(line, "\t", tokens);

        t.testID  = tokens[FTestID];
        t.geneID  = tokens[FGeneID];        
        t.fpkm_1  = stof(tokens[FFPKM_1]);
        t.fpkm_2  = stof(tokens[FFPKM_2]);
        t.status  = tok2Status.at(tokens[FStatus]);
        t.logFold = stof(tokens[FLogFold]);
        t.stats   = stof(tokens[FTestStats]);

        t.p = SS::P(stof(tokens[FPValue]));
        t.q = SS::P(stof(tokens[FQValue]));

        if (t.status != TrackingStatus::HIData)
        {
            f(t, p);
        }
    }
}


