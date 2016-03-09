#ifndef PARSER_SLEUTH_HPP
#define PARSER_SLEUTH_HPP

#include "data/dtest.hpp"
#include "data/tokens.hpp"
#include "data/standard.hpp"
#include "stats/analyzer.hpp"

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

        static void parse(const FileName &file, std::function<void (const DiffTest &, const ParserProgress &)> f)
        {
            Reader r(file);
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks;
            
            // We'll need it for checking sequins
            //const auto &ref = Standard::instance().r_trans;

            while (r.nextLine(line))
            {
                Tokens::split(line, ",", toks);

                DiffTest t;

                if (p.i)
                {
                    t.id  = toks[Field::TargetID];
                    t.cID = "chrT"; // TODO: ref.match(t.id) ? ChrT : Endo;

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
                        t.p = stod(toks[Field::PValue]);
                        
                        // Probability controlled for multi-testing
                        t.q = stod(toks[Field::QValue]);
                    }

                    f(t, p);
                }

                p.i++;
            }
        }
    };
}

#endif
