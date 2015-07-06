#ifndef GI_PARSER_FUSION_HPP
#define GI_PARSER_FUSION_HPP

#include <functional>
#include "data/locus.hpp"
#include "data/reader.hpp"
#include "data/biology.hpp"

namespace Spike
{
    struct ParserFusion
    {
        struct Fusion
        {
            inline operator Locus() const
            {
                return Locus(std::min(start_1, start_1), std::max(start_1, start_1));
            }

            // The name of the first chromosome
            std::string chr_1;
            
            // The name of the second chromosome
            std::string chr_2;

            // Position of the first chromosome where the fusion occurs
            BasePair start_1;
            
            // Position of the second chromosome where the fusion occurs
            BasePair start_2;
            
            // Orientation for the first chromosome
            Strand dir_1;
            
            // Orientation for the second chromosome
            Strand dir_2;
        };

        typedef std::function<void(const Fusion &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif