#ifndef PARSER_HTSEQ_COUNT_HPP
#define PARSER_HTSEQ_COUNT_HPP

#include "parsers/parser_csv.hpp"

namespace Anaquin
{
    struct ParserHTSeqCount
    {
        typedef unsigned Count;
        
        struct Sample
        {
            // Eg: ENSG00000000419.12
            std::string id;

            // Eg: 71
            Count count;
        };
        
        struct Samples
        {
            // Eg: "ENSG00000000419.12"
            std::string id;

            // Eg: 34, 56
            std::vector<Count> counts;
        };
        
        /*
         * Parse multiple files simultaneously. The files are assumed be sorted.
         */
        
        static void parse(const std::vector<Reader> &rs, std::function<void (const Samples &, const ParserProgress &)> f)
        {
            ParserProgress p;            
            std::vector<std::string> toks;

            Samples s;
            s.counts.resize(rs.size());

            while (true)
            {
                s.id.clear();
                
                for (auto i = 0; i < rs.size(); i++)
                {
                    if (!rs[i].nextTokens(toks, "\t"))
                    {
                        return;
                    }
                    
                    if (i && toks[0] != s.id)
                    {
                        throw std::runtime_error("Files are not sorted.");
                    }
                    
                    // Eg: ENSG00000000457.13
                    s.id = toks[0];

                    // Eg: 49
                    s.counts[i] = stoi(toks[1]);
                }
                
                f(s, p);
                p.i++;
            }
        }
        
        static void parse(const Reader &r, std::function<void (const Sample &, const ParserProgress &)>  f)
        {
            /*
             * Count table is simply a CSV file with two columns:
             *
             *      ENSG00000000419.12	71
             *      ENSG00000000457.13	49
             */
            
            Sample s;
            
            ParserCSV::parse(r, [&](const ParserCSV::Data &d, const ParserProgress &p)
            {
                if (d.size() != 2)
                {
                    throw std::runtime_error("Invalid file for HTSeqCount. Please try and check again.");
                }
                
                // Eg: ENSG00000000457.13
                s.id = d[0];
                
                // Eg: 49
                s.count = stoi(d[1]);

                f(s, p);
            }, "\t");
        }
    };
}

#endif
