#ifndef PARSER_STAR_FUSION_HPP
#define PARSER_STAR_FUSION_HPP

#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "parsers/parser.hpp"
#include "stats/analyzer.hpp"
#include "FusQuin/FUSQuin.hpp"

namespace Anaquin
{
    struct ParserStarFusion
    {
        enum Field
        {
            FusionName,
            JunctionReads,
            SpanningFrags,
            LeftGene,
            LeftBreakpoint,
            LeftDistFromRefExonSplice,
            RightGene,
            RightBreakpoint,
            RightDistFromRefExonSplice,
        };

        typedef CalledFusion Data;

        // Parse an output file from FusionStar (eg: star-fusion.fusion_candidates.txt)
        template <typename F> static void parse(const Reader &r, F f)
        {
            std::string line;
            std::vector<Tokens::Token> toks;
            
            Data data;
            ParserProgress p;
            
            while (r.nextLine(line))
            {
                if (line[0] == '#')
                {
                    continue;
                }
                
                p.i++;
                Tokens::split(line, "\t", toks);
                
                const auto leftGene  = toks[Field::LeftGene];
                const auto rightGene = toks[Field::RightGene];

                auto parseBreak = [&](const std::string &s, std::string &chr, Base &l, Strand &o)
                {
                    /*
                     * Eg: chrT:9035684:-
                     */
                    
                    std::vector<std::string> tokens;
                    Tokens::split(s, ":", tokens);
                    
                    assert(tokens.size() == 3);
                    
                    // Eg: chrT
                    chr = tokens[0];
                    
                    // Eg: 9035684
                    l = std::stoi(tokens[1]);
                    
                    // Eg: '-'
                    o = tokens[2] == "-" ? Backward : Forward;
                };

                parseBreak(toks[Field::LeftBreakpoint],  data.cID_1, data.l1, data.s1);
                parseBreak(toks[Field::RightBreakpoint], data.cID_2, data.l2, data.s2);

                // Measured abundance
                data.reads = stoi(toks[Field::JunctionReads]);
                
                f(data, p);
            }
        }
    };
}

#endif