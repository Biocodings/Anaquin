#ifndef PARSER_VARIANT_HPP
#define PARSER_VARIANT_HPP

#include "data/tokens.hpp"
#include "data/variant.hpp"
#include "data/convert.hpp"
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
            Alt,
            ReadR,
            ReadV,
            Depth,
            QualR,
            QualV,
            PValue,
        };

        typedef CalledVariant Data;

        typedef std::function<void(const Data &, const ParserProgress &)> Functor;

        static bool isVariant(const Reader &r)
        {
            std::string line;
            std::vector<Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                
                if (toks.size() == 10      &&
                    toks[0]  == "ChrID"    &&
                    toks[1]  == "Position" &&
                    toks[2]  == "Ref"      &&
                    toks[3]  == "Alt"      &&
                    toks[4]  == "ReadR"    &&
                    toks[5]  == "ReadV"    &&
                    toks[6]  == "Depth"    &&
                    toks[7]  == "QualR"    &&
                    toks[8]  == "QualV"    &&
                    toks[9]  == "PValue")
                {
                    return true;
                }
            }

            return false;
        }
        
        static void parse(const Reader &r, Functor f)
        {
            Data d;
            ParserProgress p;
            
            Line line;
            std::vector<Token> toks;

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

                d.qualR = s2d(toks[QualR]);
                d.qualV = s2d(toks[QualV]);
                d.readR = s2d(toks[ReadR]);
                d.readV = s2d(toks[ReadV]);
                d.depth = d.readR + d.readV;

                d.allF  = static_cast<Proportion>(d.readV) / (d.readR + d.readV);

                /*
                 * Eg: chrIS  631340  A  A  3976  0
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
                d.alt = toks[Alt];
                
                try
                {
                    d.p = stold(toks[PValue]);
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