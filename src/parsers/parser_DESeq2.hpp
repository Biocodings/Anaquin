#ifndef PARSER_DESEQ2_HPP
#define PARSER_DESEQ2_HPP

#include "data/dtest.hpp"
#include "data/tokens.hpp"
#include "data/standard.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserDESeq2
    {
        typedef enum
        {
            Name,
            BaseName,
            Log2Fold,
            Log2FoldSE,
            Stat,
            PValue,
            QValue
        } DESeq2Field;

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
                 * DESeq2 is an R-package, it has no output directly relevant to Anaquin. However, we'd expect
                 * the users write the resulting differential table to a text file. It's format would like this:
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
                    t.id = toks[DESeq2Field::Name];
                    
                    /*
                     * DESeq2 wouldn't give the name of the chromosome, only the name of the gene would be given.
                     * We'll simply have to consult the reference annotation whether the gene belongs to the synthetic
                     * chromosome.
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

                    try
                    {
                        // Measured log-fold change
                        t.logF = stof(toks[DESeq2Field::Log2Fold]);
                        
                        // Probability under the null hypothesis
                        t.p = stod(toks[DESeq2Field::PValue]);
                        
                        // Probability controlled for multi-testing
                        t.q = stod(toks[DESeq2Field::QValue]);
                    }
                    catch (const std::invalid_argument &)
                    {
                        continue;
                    }
                    
                    t.status = DiffTest::Status::Tested;

                    f(t, p);
                }

                p.i++;
            }
        }
    };
}

#endif
