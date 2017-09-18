#ifndef PARSER_SLEUTH_HPP
#define PARSER_SLEUTH_HPP

#include "data/dtest.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
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
        
        typedef DiffTest Data;

        static bool isSleuth(const Reader &r)
        {
            std::string line;
            std::vector<Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, ",", toks);
                
                if (toks.size() == 11             &&
                    toks[0]  == "target_id"       &&
                    toks[1]  == "pval"            &&
                    toks[2]  == "qval"            &&
                    toks[3]  == "b"               &&
                    toks[4]  == "se_b"            &&
                    toks[5]  == "mean_obs"        &&
                    toks[6]  == "var_obs"         &&
                    toks[7]  == "tech_var"        &&
                    toks[8]  == "sigma_sq"        &&
                    toks[9]  == "smooth_sigma_sq" &&
                    toks[10] == "final_sigma_sq")
                {
                    return true;
                }
            }
            
            return false;
        }
        
        template <typename F> static void parse(const Reader &r, F f)
        {
            const auto &ref = Standard::instance().r_rna;
            
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks;

            while (r.nextLine(line))
            {
                Tokens::split(line, ",", toks);
                
                Data t;
                
                if (p.i)
                {
                    t.iID = toks[Field::TargetID];
                    
                    // Can we match the isoform to sequins?
                    auto isChrIS = ref.seqsL1().count(t.iID);
                    
                    t.cID = isChrIS ? ChrIS() : "geno";
                    t.gID = ""; // TODO: isChrIS ? ref.s2g(t.iID) : "";
                    
                    if (toks[Field::PValue] == "NA" || toks[Field::QValue] == "NA")
                    {
                        t.status = DiffTest::Status::NotTested;
                    }
                    else
                    {
                        t.status = DiffTest::Status::Tested;
                        
                        t.mean = stod(toks[Field::MeanObs]);
                        
                        // Measured log-fold change
                        t.logF_ = stod(toks[Field::B]);
                        
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
