#ifndef PARSER_VARSCAN_HPP
#define PARSER_VARSCAN_HPP

#include <iostream>
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
            PValue,
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

        static bool isVarScan(const Reader &r)
        {
            std::string line;
            std::vector<Tokens::Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                
                if (toks.size() == 19         &&
                    toks[0]  == "Chrom"       &&
                    toks[1]  == "Position"    &&
                    toks[2]  == "Ref"         &&
                    toks[3]  == "Cons"        &&
                    toks[4]  == "Reads1"      &&
                    toks[5]  == "Reads2"      &&
                    toks[6]  == "VarFreq"     &&
                    toks[7]  == "Strands1"    &&
                    toks[8]  == "Strands2"    &&
                    toks[9]  == "Qual1"       &&
                    toks[10] == "Qual2"       &&
                    toks[11] == "Pvalue"      &&
                    toks[12] == "MapQual1"    &&
                    toks[13] == "MapQual2"    &&
                    toks[14] == "Reads1Plus"  &&
                    toks[15] == "Reads1Minus" &&
                    toks[16] == "Reads2Plus"  &&
                    toks[17] == "Reads2Minus" &&
                    toks[18] == "VarAllele")
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

                d.allF  = stod(toks[VarFreq]);
                
                /*
                 * TODO: VarScan does give quality but it's not the same as quality in VCF.
                 *       VarScan gives quality score for the reference and alternative allele.
                 */
                
                const auto readR = stod(toks[Reads1]);
                const auto readV = stod(toks[Reads2]);
                
                d.readR = readR;
                d.readV = readV;
                d.depth = d.readR + d.readV;
                
                if (!readV)
                {
                    continue;
                }
                
                /*
                 * Is this an insertion?
                 */
                
                const auto isInsert = toks[Cons].find("*/+") != std::string::npos;
                const auto isDelete = toks[Cons].find("*/-") != std::string::npos;

                auto clean = [&](const std::string &s)
                {
                    auto x = s;
                    
                    boost::replace_all(x, "*", "");
                    boost::replace_all(x, "/", "");
                    boost::replace_all(x, "-", "");
                    boost::replace_all(x, "+", "");
                    
                    return x;
                };
                
                if (isInsert)
                {
                    /*
                     * Eg: A * /+CT	+CT
                     */
                    
                    d.ref = clean(toks[Ref]);
                    d.alt = clean(toks[Ref] + toks[VarAllele]);
                }
                else if (isDelete)
                {
                    /*
                     * Eg: G * /-CTTCCTCTTTC CTTCCTCTTTC
                     */
                    
                    d.ref = clean(toks[Ref] + toks[VarAllele]);
                    d.alt = clean(toks[Ref]);
                }
                else
                {
                    d.ref = clean(toks[Ref]);
                    d.alt = clean(toks[VarAllele]);
                }
                
                assert(d.ref.find('+') == std::string::npos);
                assert(d.alt.find('+') == std::string::npos);
                
                try
                {
                    d.p = stold(toks[PValue]);
                }
                catch (...)
                {
                    d.p = 0.0;
                }
                
                assert(d.p >= 0 && d.p <= 1.0);

                f(d, p);
            }
        }
    };
}

#endif