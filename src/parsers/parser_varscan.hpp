#ifndef PARSER_VARSCAN_HPP
#define PARSER_VARSCAN_HPP

#include "data/tokens.hpp"
#include "data/variant.hpp"
#include "parsers/parser.hpp"
#include <boost/algorithm/string.hpp>

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

                d.readR = stod(toks[Reads1]);
                d.readV = stod(toks[Reads2]);
                
                // Why it's shown as a variant if no variant supporting reads?
                if (!d.readV)
                {
                    continue;
                }
                
                if (toks[Cons][0] == '*')
                {
                    /*
                     * Eg:
                     *
                     *     chrT	8288872	C	* /-CCTG	2443	367	13.05%	2	2	29	19	1.8024382198605343E-112	1	1	1228	1215	186	181	-CCTG
                     *
                     */

                    // Eg: C*/-CCTG
                    d.ref = toks[Ref] + toks[Cons];

                    // Eg: C/-CCTG
                    boost::replace_all(d.ref, "*", "");
                    
                    // Eg: C-CCTG
                    boost::replace_all(d.ref, "/", "");
                    
                    // Eg: CCCTG
                    boost::replace_all(d.ref, "-", "");
                    
                    // Eg: -CCTG
                    auto tmp = toks[VarAllele];
                    
                    // Eg: CCTG
                    boost::replace_all(tmp, "-", "");
                    
                    // Eg: C
                    boost::replace_all(d.alt = d.ref, tmp, "");
                }
                else
                {
                    d.ref = toks[Ref];
                    d.alt = toks[VarAllele];
                }
                
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