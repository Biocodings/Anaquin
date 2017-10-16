#ifndef PARSER_CDIFF_HPP
#define PARSER_CDIFF_HPP

#include "data/dtest.hpp"
#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "stats/analyzer.hpp"
#include <boost/algorithm/string/predicate.hpp>

namespace Anaquin
{
    struct ParserCDiff
    {
        enum Field
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
        
        struct Data : public DiffTest
        {
            // Test statistics
            double stats;
        };
        
        template <typename F> static void parse(const Reader &r, F f)
        {
            typedef std::string TrackID;

            static const std::map<TrackID, DiffTest::Status> tok2Status =
            {
                { "OK",     DiffTest::Status::Tested    },
                { "HIDATA", DiffTest::Status::NotTested },
                { "NOTEST", DiffTest::Status::NotTested },
                { "FAIL"  , DiffTest::Status::NotTested },
            };
            
            Data t;
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks, temp;
            
            while (r.nextLine(line))
            {
                p.i++;
                
                if (boost::starts_with(line, "test_id"))
                {
                    continue;
                }
                
                Tokens::split(line, "\t", toks);
                
                t.gID    = toks[FGeneID];
                t.iID    = toks[FTestID];
                t.status = tok2Status.at(toks[FStatus]);
                t.logF_  = stof(toks[FLogFold]);
                t.stats  = stof(toks[FTestStats]);
                
                // Eg: chrIS:1082119-1190836
                Tokens::split(toks[FLocus], ":", temp);
                
                // Eg: chrIS
                t.cID = temp[0];
                
                t.p = stold(toks[FPValue]);
                t.q = stold(toks[FQValue]);
                
                if (t.status != DiffTest::Status::NotTested)
                {
                    f(t, p);
                }
            }
        }
        
        static bool isCDiff(const Reader &r, unsigned &nTrans, unsigned &nGenes);
    };
}

#endif
