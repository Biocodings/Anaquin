#ifndef PARSER_EXP_HPP
#define PARSER_EXP_HPP

#include <map>
#include <cctype>
#include <vector>
#include "data/data.hpp"
#include "data/reader.hpp"

namespace Anaquin
{
    struct ParserExp
    {
        struct Sample
        {
            // File name for the first paired-end
            FileName first;
            
            // File anme for the second paired-end
            FileName second;
        };
        
        struct Experiment
        {
            std::map<Mixture, std::vector<Sample>> samps;
        };

        static Experiment parse(const Reader &r) throw()
        {
            Experiment exp;
            
            protectParse("Experiment format", [&]()
            {
                std::vector<std::string> toks;

                for (auto i = 0; r.nextTokens(toks, "\t"); i++)
                {
                    if (toks.size() != 3)
                    {
                        throw BadFormatException("Invalid line: " + r.lastLine());
                    }
                    else if (!i)
                    {
                        continue;
                    }
                    
                    std::transform(toks[0].begin(), toks[0].end(), toks[0].begin(), ::toupper);
                    
                    if (toks[0] != "A" && toks[0] != "B")
                    {
                        throw BadFormatException("Invalid mixture: " + toks[0]);
                    }
                    
                    Sample samp;
                    
                    samp.first  = toks[1];
                    samp.second = toks[2];

                    if (toks[0] == "A") { exp.samps[Mixture::Mix_1].push_back(samp); }
                    else                { exp.samps[Mixture::Mix_2].push_back(samp); }
                }
            });

            return exp;
        }
    };
}

#endif