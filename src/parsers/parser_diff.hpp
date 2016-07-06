/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef PARSER_DIFF_HPP
#define PARSER_DIFF_HPP

#include "data/data.hpp"
#include "data/tokens.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserDiff
    {
        typedef enum
        {
            ChrID,
            GeneID,
            IsoformID,
            Log2Fold,
            Log2FoldSE,
            PValue,
            QValue,
            Mean
        } Field;

        struct Data
        {
            ::Anaquin::ChrID cID;

            ::Anaquin::GeneID gID;
            ::Anaquin::IsoformID iID;
            
            double logF;
            
            // Not always reported (eg: Cuffdiff)
            double logFSE = NAN;

            Probability p;
            
            // Adjusted p-values
            Probability q;
            
            // Eg: Basemean in DESeq2
            double mean;

            inline bool hasSE() const
            {
                return !isnan(logFSE);
            }
        };

        template <typename F> static void parse(const FileName &file, F f)
        {
            Reader r(file);
            ParserProgress p;
            std::vector<Tokens::Token> toks;

            std::string line;
            
            while (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                Data x;

                if (p.i)
                {
                    x.cID    = toks[Field::ChrID];
                    x.gID    = toks[Field::GeneID];
                    x.iID    = toks[Field::IsoformID] != "-" ? toks[Field::GeneID] : "";
                    x.p      = s2p(toks[Field::PValue]);
                    x.q      = s2p(toks[Field::QValue]);
                    x.mean   = s2d(toks[Field::Mean]);
                    x.logF   = s2d(toks[Field::Log2Fold]);
                    x.logFSE = s2d(toks[Field::Log2FoldSE]);

                    toString("");
                    f(x, p);
                }

                p.i++;
            }
        }
        
        static bool isIsoform(const Reader &r)
        {
            std::vector<Tokens::Token> toks;
            
            if (r.nextTokens(toks, "\t"))
            {
                return toks[1] == "IsoID";
            }
            
            return true;
        }
    };
}

#endif
