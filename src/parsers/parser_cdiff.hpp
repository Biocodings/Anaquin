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
        typedef std::string TrackID;
        
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
        
        struct Data : public DiffTest
        {
            TrackID tID;
            
            // Test statistics
            double stats;
        };

        template <typename F> static void parse(const FileName &file, F f)
        {
            static const std::map<ParserCDiff::TrackID, DiffTest::Status> tok2Status =
            {
                { "OK",     DiffTest::Status::Tested    },
                { "HIDATA", DiffTest::Status::NotTested },
                { "NOTEST", DiffTest::Status::NotTested },
            };

            Reader i(file);
            
            Data t;
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
                t.tID    = toks[FTestID];
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
    };    
}

#endif
