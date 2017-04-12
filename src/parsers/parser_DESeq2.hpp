#ifndef PARSER_DESEQ2_HPP
#define PARSER_DESEQ2_HPP

#include "data/dtest.hpp"
#include "data/tokens.hpp"
#include "tools/tools.hpp"
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
        } Field;
        
        typedef DiffTest Data;

        static bool isDESeq2(const Reader &r)
        {
            std::string line;
            std::vector<Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, ",", toks);
                
                if (toks.size() == 7             &&
                    toks[1]  == "baseMean"       &&
                    toks[2]  == "log2FoldChange" &&
                    toks[3]  == "lfcSE"          &&
                    toks[4]  == "stat"           &&
                    toks[5]  == "pvalue"         &&
                    toks[6]  == "padj")
                {
                    return true;
                }
            }
            
            return false;
        }

        template <typename F> static void parse(const FileName &file, F f)
        {
            Reader r(file);
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks;
            
            // We'll need it for checking sequin genes
            const auto &s = Standard::instance().r_rna;
            
            while (r.nextLine(line))
            {
                Tokens::split(line, ",", toks);

                Data t;

                if (p.i)
                {
                    t.gID = toks[Field::Name];
                    
                    if (t.gID == "R2_33")
                    {
                        t.gID = t.gID;
                    }
                    
                    
                    /*
                     * DESeq2 wouldn't give the chromosomes, only the gene names.
                     * We have to consult the reference annotation to make a decision.
                     */
                    
                    t.cID = s.findGene(ChrIS, t.gID) ? ChrIS : Geno;
                    
                    /*
                     * Eg: ENSG00000000003.14,0,NA,NA,NA,NA,NA
                     */
                    
                    
                    if (toks[Field::PValue] == "NA" || toks[Field::Log2Fold] == "NA")
                    {
                        t.status = DiffTest::Status::NotTested;
                    }
                    else
                    {
                        t.status = DiffTest::Status::Tested;

                        // Normalized average counts
                        t.mean = s2d(toks[Field::BaseMean]);
                        
                        // Measured log-fold change
                        t.logF_ = s2d(toks[Field::Log2Fold]);

                        t.samp1 = NAN;
                        t.samp2 = NAN;

                        // Standard error for the log-fold change
                        t.logFSE = s2d(toks[Field::Log2FoldSE]);
                        
                        t.p = stold(toks[Field::PValue]);
                        
                        try
                        {
                            // Not always available, but we can still proceed if we have p-value
                            t.q = stold(toks[Field::QValue]);
                        }
                        catch (...)
                        {
                            t.q = NAN;
                        }
                    }

                    f(t, p);
                }

                p.i++;
            }
        }
    };
}

#endif
