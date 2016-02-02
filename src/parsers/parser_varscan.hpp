#ifndef PARSER_VARSCAN_HPP
#define PARSER_VARSCAN_HPP

#include "data/tokens.hpp"
#include "data/variant.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserVarScan
    {
        enum Field
        {
            Chrom,
            Position,
            Ref,
            Cons,
            Reads1,
            Reads2,
            VarFreq,
            Strands1,
            Strands2,
            Qual1,
            Qual2,
            Pvalue,
            MapQual1,
            MapQual2,
            Reads1Plus,
            Reads1Minus,
            Reads2Plus,
            Reads2Minus,
            VarAllele
        };

        typedef CalledVariant Data;

        typedef std::function<void(const Data &, const ParserProgress &)> Functor;

        static void parse(const Reader &r, Functor f)
        {
            Data d;
            ParserProgress p;
            
            Line line;
            std::vector<Tokens::Token> toks;

            while (r.nextLine(line))
            {
                if (p.i++ == 0)
                {
                    continue;
                }
                
                Tokens::split(line, "\t", toks);
                
                d.chrID = toks[Chrom];

                // Always start and end at the same position
                d.l = Locus(stod(toks[Position]), stod(toks[Position]));

                d.ref = toks[Ref];
                d.alt = toks[VarAllele];
                
                d.readR = stod(toks[Reads1]);
                d.readV = stod(toks[Reads2]);
                
                try
                {
                    d.pval = stod(toks[Pvalue]);
                }
                catch (...)
                {
                    d.pval = 0.0;
                }

                f(d, p);
            }
        }
    };
}

#endif