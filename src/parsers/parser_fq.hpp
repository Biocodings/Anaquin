#ifndef PARSER_FQ_HPP
#define PARSER_FQ_HPP

#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserFQ
    {
        struct Data
        {
            std::string name, info, seq, opt, qual;
        };

        template <typename F> static void parse(const Reader &r, F f)
        {
            Data x;
            ParserProgress p;
            std::vector<std::string> toks;
            
            while (r.nextLine(x.name))
            {
                if (!r.nextLine(x.seq))  { throw "Error. Is this a valid FASTQ file?"; }
                if (!r.nextLine(x.opt))  { throw "Error. Is this a valid FASTQ file?"; }
                if (!r.nextLine(x.qual)) { throw "Error. Is this a valid FASTQ file?"; }

                Tokens::split(x.name, " ", toks);
                x.name = toks[0];

                // TODO:...
                if (toks.size() > 1)
                {
                    x.info = toks[1];
                }
                
                f(x, p);
                p.i++;
            }
        }
    };
}

#endif