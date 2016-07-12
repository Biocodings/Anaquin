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
            LogFC,
            LogCPM,
            PValue
        } Field;

        static bool isEdgeR(const Reader &r)
        {
            std::string line;
            std::vector<Tokens::Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);

                if (toks.size() == 3     &&
                    toks[0]  == "logFC"  &&
                    toks[1]  == "logCPM" &&
                    toks[2]  == "PValue")
                {
                    return true;
                }
            }
            
            return false;
        }
        
        typedef DiffTest Data;

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
                    t.gID = toks[Field::Name];

                    /*
                     * edgeR wouldn't give the chromosome name, only the name of the gene would be given.
                     * We have to consult the reference annotation to make a decision.
                     */
                    
                    if (s.findGene(ChrT, t.gID))
                    {
                        t.cID = ChrT;
                    }
                    else
                    {
                        t.cID = Geno;
                    }
                    
                    if (toks[Field::PValue] == "NA" || toks[Field::LogFC] == "NA")
                    {
                        t.status = DiffTest::Status::NotTested;
                    }
                    else
                    {
                        t.status = DiffTest::Status::Tested;
                        
                        // Normalized average counts
                        //t.baseMean = stod(toks[Field::BaseMean]);
                        
                        // Measured log-fold change
                        t.logF = stod(toks[Field::LogFC]);
                        
                        // Standard error for the log-fold change
                        //t.logFSE = stod(toks[Field::Log2FoldSE]);
                        
                        // Probability under the null hypothesis
                        t.p = stold(toks[Field::PValue]);
                        
                        // Probability adjusted for multi-testing
                        //t.q = stold(toks[Field::QValue]);
                    }
                    
                    f(t, p);
                }
                
                p.i++;
            }
        }        
    };
}

#endif
