#ifndef GI_PARSER_TOP_FUSION_HPP
#define GI_PARSER_TOP_FUSION_HPP

#include <functional>
#include "data/locus.hpp"
#include "data/reader.hpp"
#include "data/biology.hpp"

namespace Anaquin
{
    struct ParserTopFusion
    {
        struct Fusion
        {
            inline operator Locus() const
            {
                return Locus(l1, l2);
            }

            // The name of the first chromosome
            std::string chr_1;
            
            // The name of the second chromosome
            std::string chr_2;

            // Position of the first chromosome where the fusion occurs
            Base l1;
            
            // Position of the second chromosome where the fusion occurs
            Base l2;
            
            // Orientation for the first chromosome
            Strand s1;
            
            // Orientation for the second chromosome
            Strand s2;

            // Number of reads that span the fusion
            Reads reads;
        };

        typedef std::function<void(const Fusion &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif