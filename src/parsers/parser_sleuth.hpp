#ifndef PARSER_SLEUTH_HPP
#define PARSER_SLEUTH_HPP

#include "data/dtest.hpp"
#include "data/tokens.hpp"
#include "data/standard.hpp"
#include "stats/analyzer.hpp"
#include "tools/gtf_data.hpp"

namespace Anaquin
{
    struct ParserSleuth
    {
        typedef enum
        {
            TargetID,
            PValue,
            QValue,
            B,
            SE_B,
            MeanObs,
            VarObs,
            TechVar,
            SigmaSQ,
            SmoothSigmaSQ,
            FinalSigmaSQ,
        } Field;

        typedef DiffTest Data;
        
        static void parse(const FileName &file, std::function<void (const Data &, const ParserProgress &)> f)
        {
            const auto &r = Standard::instance().r_trans;

            Reader read(file);
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks;
            
            while (read.nextLine(line))
            {
                Tokens::split(line, ",", toks);

                Data t;

                if (p.i)
                {
                    t.iID = toks[Field::TargetID];
                    
                    // Can we match to the sequins?
                    auto m = r.match(t.iID);

                    t.cID = m ? ChrT : Geno;
                    t.gID = m ? r.s2g(t.iID) : "";
                    
                    if (toks[Field::PValue] == "NA" || toks[Field::QValue] == "NA")
                    {
                        t.status = DiffTest::Status::NotTested;
                    }
                    else
                    {
                        t.status = DiffTest::Status::Tested;

                        t.baseMean = stod(toks[Field::MeanObs]);

                        // Measured log-fold change
                        t.logF = stod(toks[Field::B]);

                        // Standard error for the log-fold change
                        t.logFSE = stod(toks[Field::SE_B]);
                        
                        // Probability under the null hypothesis
                        t.p = stold(toks[Field::PValue]);
                        
                        // Probability controlled for multi-testing
                        t.q = stold(toks[Field::QValue]);
                    }

                    f(t, p);
                }

                p.i++;
            }
        }
    };
}

#endif
