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
            BaseMean,
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

                if (p.i)
                {
                    t.id = toks[DESeq2Field::Name];
                    
                    /*
                     * DESeq2 wouldn't give the chromosome name, only the name of the gene would be given.
                     * We have to consult the reference annotation to make a decision.
                     */
                    
                    if (s.findGene(ChrT, t.id))
                    {
                        t.cID = ChrT;
                    }
                    else
                    {
                        t.cID = Geno;
                    }
                    
                    /*
                     * Handle the NA case:
                     *
                     *     ENSG00000000003.14,0,NA,NA,NA,NA,NA
                     */
                    
                    if (toks[DESeq2Field::PValue]   == "NA" ||
                        toks[DESeq2Field::QValue]   == "NA" ||
                        toks[DESeq2Field::Log2Fold] == "NA")
                    {
                        t.status = DiffTest::Status::NotTested;
                    }
                    else
                    {
                        t.status = DiffTest::Status::Tested;

                        // Normalized average counts
                        t.baseMean = stod(toks[DESeq2Field::BaseMean]);
                        
                        // Measured log-fold change
                        t.logF = stod(toks[DESeq2Field::Log2Fold]);

                        // Standard error for the log-fold change
                        t.logFSE = stod(toks[DESeq2Field::Log2FoldSE]);
                        
                        // Probability under the null hypothesis
                        t.p = stold(toks[DESeq2Field::PValue]);
                        
                        // Probability adjusted for multi-testing
                        t.q = stold(toks[DESeq2Field::QValue]);                        
                    }

                    f(t, p);
                }

                p.i++;
            }
        }
    };
}

#endif
