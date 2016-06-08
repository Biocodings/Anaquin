#ifndef PARSER_VARIANT_HPP
#define PARSER_VARIANT_HPP

#include <iostream>
#include "data/tokens.hpp"
#include "data/variant.hpp"
#include "parsers/parser.hpp"
#include <boost/algorithm/string.hpp>

namespace Anaquin
{
    struct ParserVariant
    {
        enum Field
        {
            Chrom,
            Position,
            Ref,
            Allele,
            Reads1,
            Reads2,
            Qual1,
            Qual2,
            Pvalue,
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

                d.qualR = stod(toks[Qual1]);
                d.qualV = stod(toks[Qual2]);
                d.readR = stod(toks[Reads1]);
                d.readV = stod(toks[Reads2]);
                d.allF  = static_cast<Proportion>(d.readV) / (d.readR + d.readV);

                /*
                 * Eg:
                 *
                 *   chrT  631340  A  A  3976  0
                 *
                 * In the example, the sixth column is the reads for the variant allele.
                 *
                 * Why it's shown as a variant if no variant supporting the variant allele?
                 */
                
                if (!d.readV)
                {
                    continue;
                }
                
                d.ref = toks[Ref];
                d.alt = toks[Allele];
                
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