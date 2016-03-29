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
            Qual2,       // 10
            Pvalue,
            MapQual1,
            MapQual2,
            Reads1Plus,
            Reads1Minus,
            Reads2Plus,
            Reads2Minus,
            VarAllele    // 18
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
                
                d.cID = toks[Chrom];

                // Always start and end at the same position
                d.l = Locus(stod(toks[Position]), stod(toks[Position]));

                d.ref = toks[Ref];

                d.readR = stod(toks[Reads1]);
                d.readV = stod(toks[Reads2]);

                // Why it's shown as a variant if there's no variant supporting reads?
                if (!d.readV)
                {
                    continue;
                }

                d.alt = toks[VarAllele];
                
                try
                {
                    d.p = stod(toks[Pvalue]);
                }
                catch (...)
                {
                    d.p = 0.0;
                }

                f(d, p);
            }
        }
    };
}

#endif