/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

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
            // Test statistics
            double stats;
        };

        static bool isTracking(const Reader &r)
        {
            std::string line;
            std::vector<Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                
                if (toks.size() == 14               &&
                    toks[0]  == "test_id"           &&
                    toks[1]  == "gene_id"           &&
                    toks[2]  == "gene"              &&
                    toks[3]  == "locus"             &&
                    toks[4]  == "sample_1"          &&
                    toks[5]  == "sample_2"          &&
                    toks[6]  == "status"            &&
                    toks[7]  == "value_1"           &&
                    toks[8]  == "value_2"           &&
                    toks[9]  == "log2(fold_change)" &&
                    toks[10] == "test_stat"         &&
                    toks[11] == "p_value"           &&
                    toks[12] == "q_value"           &&
                    toks[13] == "significant")
                {
                    return true;
                }
            }

            return false;
        }
        
        template <typename F> static void parse(const FileName &file, F f)
        {
            static const std::map<ParserCDiff::TrackID, DiffTest::Status> tok2Status =
            {
                { "OK",     DiffTest::Status::Tested    },
                { "HIDATA", DiffTest::Status::NotTested },
                { "NOTEST", DiffTest::Status::NotTested },
                { "FAIL"  , DiffTest::Status::NotTested },
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

                t.gID    = toks[FGeneID];
                t.iID    = toks[FTestID];
                t.status = tok2Status.at(toks[FStatus]);
                t.logF_   = stof(toks[FLogFold]);
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
    };    
}

#endif
