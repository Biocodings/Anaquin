#ifndef PARSER_EDGER_HPP
#define PARSER_EDGER_HPP

#include "data/dtest.hpp"
#include "data/tokens.hpp"

namespace Anaquin
{
    struct ParserEdgeR
    {
        typedef enum
        {
            Name,
            BaseMean,
            Log2Fold,
            Log2FoldSE,
            Stat,
            PValue,
            QValue
        } Field;
        
        static void parse(const FileName &file, std::function<void (const DiffTest &, const ParserProgress &)> f)
        {
            Reader r(file);
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks;
            
            // We'll need it for checking sequin genes
            const auto &s = Standard::instance().r_trans;
            
            while (r.nextLine(line))
            {
                Tokens::split(line, ",", toks);
                
                DiffTest t;
                
                /*
                 * EdgeR is an R-package, it has no output directly relevant to Anaquin. We'd expect the user
                 * write the resulting differential table to a text file.
                 *
                 *     Column 1: Name
                 *     Column 2: Base mean
                 *     Column 3: Log2 Fold
                 *     Column 4: Log2 Fold SE
                 *     Column 5: Stat
                 *     Column 6: P-value
                 *     Column 7: Q-value
                 */
                
                if (p.i)
                {
                    t.id = toks[Field::Name];
                    
                    /*
                     * EdgeR wouldn't give the name of the chromosome, only the name of the gene would be given.
                     * We'll have to consult the reference annotation whether the gene belongs to the synthetic.
                     */
                    
                    if (s.findGene(ChrT, t.id))
                    {
                        t.cID = ChrT;
                    }
                    
                    /*
                     * In differential analysis, it doens't matter exactly where the endogenous genes come from. They'll
                     * all have the same label.
                     */
                    
                    else
                    {
                        t.cID = Endo;
                    }
                    
                    /*
                     * Handle the NA case:
                     *
                     *     ENSG00000000003.14,0,NA,NA,NA,NA,NA
                     */
                    
                    if (toks[Field::PValue]   == "NA" ||
                        toks[Field::QValue]   == "NA" ||
                        toks[Field::Log2Fold] == "NA")
                    {
                        t.status = DiffTest::Status::NotTested;
                    }
                    else
                    {
                        t.status = DiffTest::Status::Tested;
                        
                        // Normalized average counts
                        t.baseMean = stod(toks[Field::BaseMean]);
                        
                        // Measured log-fold change
                        t.logF = stod(toks[Field::Log2Fold]);
                        
                        // Standard error for the log-fold change
                        t.logFSE = stod(toks[Field::Log2FoldSE]);
                        
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
