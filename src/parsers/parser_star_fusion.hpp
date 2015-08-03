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
            inline operator Locus() const
            {
                return Locus(std::min(start_1, start_1), std::max(start_1, start_1));
            }
            
            // Chromosome for the left and right
            ChromoID chr_1, chr_2;

            // Strand for the left and right
            Strand strand_1, strand_2;

            // Where the fusion occurs
            BasePair start_1, start_2;
            
            // Number of reads that span the fusion
            Reads reads;
        };

        typedef std::function<void (const Fusion &, const ParserProgress &)> Functor;

        // Parse an output file from FusionStar
        static void parse(const Reader &, Functor);
    };
}

#endif