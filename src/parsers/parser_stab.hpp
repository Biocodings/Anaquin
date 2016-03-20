#ifndef PARSER_STAB_HPP
#define PARSER_STAB_HPP

#include "data/tokens.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserSTab
    {
        enum VCFField
        {
            Chromo,
            Start,
            End,
            Strand,
            Motif,
            Annotation,
            UniqueReads,
            MultiReads,
            Overhang,
        };

        struct Chimeric
        {
            ChrID id;
            
            // Base of the intron
            Locus l;
            
            // Number of uniquely mapping reads
            Counts unique;
            
            // Number of multi-mapping reads
            Counts multi;
        };
        
        // Eg: SJ.out.tab
        template <typename F> static void parse(const Reader &r, F f)
        {
            /*
             * 1: Chromosome
             * 2: First base of the intron (1-based)
             * 3: Last base of the intron (1-based)
             * 4: Strand (0: undefined, 1: +, 2: -)
             * 5: Intron motif
             * 6: Annotation
             * 7: Number of uniquely mapping reads crossing the junction
             * 8: Number of multi-mapping reads crossing the junction
             * 9: Maximum spliced alignment overhang
             */

            std::string line;
            std::vector<std::string> fields;

            ParserProgress p;
            
            while (r.nextLine(line))
            {
                p.i++;
                
                Tokens::split(line, "\t", fields);
                Chimeric c;

                c.id = fields[Chromo];
                c.l  = Locus(stod(fields[Start]), stod(fields[End]));

                c.unique = stod(fields[UniqueReads]);
                c.multi  = stod(fields[MultiReads]);

                f(c, p);
            }
        }
    };
}

#endif