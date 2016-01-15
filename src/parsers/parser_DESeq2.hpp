#ifndef PARSER_DESEQ2_HPP
#define PARSER_DESEQ2_HPP

#include "data/dtest.hpp"
#include "data/tokens.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserDESeq2
    {
        typedef enum
        {
            Name    = 0,
            PVal    = 2,
            QVal    = 2,
            LogFold = 3,
        } DESeq2Field;

        static void parse(const FileName &file, std::function<void (const DiffTest &, const ParserProgress &)> f)
        {
            Reader r(file);
            ParserProgress p;
            
            std::string line;
            std::vector<std::string> toks;
            
            while (r.nextLine(line))
            {
                Tokens::split(line, ",", toks);
                
                DiffTest t;
                
                /*
                 * DESeq2 is an R-package, it has no output directly relevant to Anaquin. However, we'd expect
                 * the users write the resulting differential table to a text file. It's format would like this:
                 *
                 *     Column 1: Gene name
                 *     Column 3: Q-value
                 */
                
                if (p.i)
                {
                    t.cID = ChrT; // TODO: Fix this
                    
                    // Eg: R1_43
                    t.id = toks[DESeq2Field::Name];

                    // Measured log-fold change
                    t.logF = stof(toks[DESeq2Field::LogFold]);

                    t.p = stod(toks[DESeq2Field::PVal]);
                    t.q = stod(toks[DESeq2Field::QVal]);

                    t.status = DiffTest::Status::Tested;
                    
                    f(t, p);
                }

                p.i++;
            }
        }
    };
}

#endif
