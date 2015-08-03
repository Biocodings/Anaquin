#ifndef GI_PARSER_STAR_FUSION_HPP
#define GI_PARSER_STAR_FUSION_HPP

#include <functional>
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserStarFusion
    {
        struct Fusion
        {
            // Chromosome for the left and right
            ChromoID l_chr, r_chr;

            // Strand for the left and right
            Strand l_strand, r_strand;

            // Where the fusion occurs
            BasePair l_break, r_break;
        };

        typedef std::function<void (const Fusion &, const ParserProgress &)> Functor;

        static void parse(const Reader &, Functor);
    };
}

#endif
