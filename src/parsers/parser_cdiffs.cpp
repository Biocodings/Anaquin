#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "parsers/parser_cdiffs.hpp"
#include <boost/algorithm/string/predicate.hpp>

using namespace Anaquin;

enum TrackingField
{
    FTestID    = 0,
    FGeneID    = 1,
    FLocus     = 3,
    FStatus    = 6,
    FFPKM_1    = 7,
    FFPKM_2    = 8,
    FLogFold   = 9,
    FTestStats = 10,
    FPValue    = 11,
    FQValue    = 12,
};

static const std::map<TrackID, DiffTest::Status> tok2Status =
{
    { "OK",     DiffTest::Status::Tested    },
    { "HIDATA", DiffTest::Status::NotTested },
    { "NOTEST", DiffTest::Status::NotTested },
};

void ParserCDiffs::parse(const FileName &file, std::function<void (const TrackingDiffs &, const ParserProgress &)> f)
{
    Reader i(file);
    TrackingDiffs t;
    ParserProgress p;

    std::string line;
    std::vector<std::string> toks, temp;
    
    while (i.nextLine(line))
    {
        p.i++;

        if (boost::starts_with(line, "test_id"))
        {
            continue;
        }
        
        Tokens::split(line, "\t", toks);

        t.id     = toks[FGeneID];
        t.testID = toks[FTestID];
        t.fpkm_1 = stof(toks[FFPKM_1]);
        t.fpkm_2 = stof(toks[FFPKM_2]);
        t.status = tok2Status.at(toks[FStatus]);
        t.logF   = stof(toks[FLogFold]);
        t.stats  = stof(toks[FTestStats]);

        // Eg: chrT:1082119-1190836
        Tokens::split(toks[FLocus], ":", temp);
        
        // Eg: chrT
        t.cID = temp[0];
        
        t.p = stof(toks[FPValue]);
        t.q = stof(toks[FQValue]);

        if (t.status != DiffTest::Status::NotTested)
        {
            f(t, p);
        }
    }
}