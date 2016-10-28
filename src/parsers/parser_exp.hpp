#ifndef PARSER_EXP_HPP
#define PARSER_EXP_HPP

#include <map>
#include <cctype>
#include <vector>
#include "data/data.hpp"
#include "data/reader.hpp"
#include "tools/errors.hpp"

namespace Anaquin
{
    struct ParserExp
    {
        struct Sample
        {
            // File name for the first paired-end
            FileName p1;
            
            // File anme for the second paired-end
            FileName p2;
        };
        
        struct Experiment
        {
            std::map<Mixture, std::vector<Sample>> samps;
        };

        static Experiment parse(const Reader &r) throw(InvalidFormatException)
        {
            Experiment exp;

            std::vector<std::string> toks;
            
            for (auto i = 0; r.nextTokens(toks, "\t"); i++)
            {
                if (toks.size() != 3)
                {
                    throw InvalidFormatException("Invalid line: " + r.lastLine());
                }
                else if (!i)
                {
                    continue;
                }
                
                std::transform(toks[0].begin(), toks[0].end(), toks[0].begin(), ::toupper);
                
                if (toks[0] != "A" && toks[0] != "B")
                {
                    throw InvalidFormatException("Invalid mixture: " + toks[0]);
                }
                
                Sample samp;
                
                samp.p1 = toks[1];
                samp.p2 = toks[2];
                
                if (toks[0] == "A") { exp.samps[Mixture::Mix_1].push_back(samp); }
                else                { exp.samps[Mixture::Mix_2].push_back(samp); }
            }

            return exp;
        }
    };
}

#endif
