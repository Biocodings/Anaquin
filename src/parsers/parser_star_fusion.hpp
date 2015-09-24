#ifndef PARSER_STAR_FUSION_HPP
#define PARSER_STAR_FUSION_HPP

#include <functional>
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserStarFusion
    {
        struct Fusion
        {
            inline operator Locus() const
            {
                return Locus(l1, l2);
            }

            // Chromosome for the left and right
            ChromoID chr_1, chr_2;

            // Strand for the left and right
            Strand s1, s2;

            // Where the fusion occurs
            Base l1, l2;
            
            // Number of reads that span the fusion
            Reads reads;
        };

        typedef std::function<void (const Fusion &, const ParserProgress &)> Functor;

        // Parse an output file from FusionStar
        static void parse(const Reader &, Functor);
    };
}

#endif