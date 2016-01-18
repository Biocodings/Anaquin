#ifndef PARSER_HTSEQ_COUNT_HPP
#define PARSER_HTSEQ_COUNT_HPP

#include "parsers/parser_csv.hpp"

namespace Anaquin
{
    struct ParserHTSeqCount
    {
        struct CountRow
        {
            // Eg: ENSG00000000419.12
            std::string id;

            // Eg: 71
            unsigned count;
        };
        
        typedef std::function<void (const CountRow &, const ParserProgress &)> Functor;

        static void parse(const Reader &r, Functor f)
        {
            /*
             * Count table is simply a CSV file with two columns:
             *
             *      ENSG00000000419.12	71
             *      ENSG00000000457.13	49
             */
            
            CountRow c;
            
            ParserCSV::parse(r, [&](const ParserCSV::Fields &fields, const ParserProgress &p)
            {
                if (fields.size() != 2)
                {
                    throw std::runtime_error("Invalid file for HTSeqCount. Please try and check again.");
                }
                
                // Eg: ENSG00000000457.13
                c.id = fields[0];
                
                // Eg: 49
                c.count = stoi(fields[1]);

                f(c, p);
            }, "\t");
        }
    };
}

#endif
