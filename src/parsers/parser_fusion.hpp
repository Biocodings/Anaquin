#ifndef GI_PARSER_FUSION_HPP
#define GI_PARSER_FUSION_HPP

#include <functional>
#include "data/types.hpp"
#include "data/reader.hpp"

namespace Spike
{
    struct ParserFusion
    {
        enum Orientation
        {
            Forward,
            Backward,
        };
        
        struct Fusion
        {
            // The name of the first chromosome
            std::string chr_1;
            
            // The name of the second chromosome
            std::string chr_2;

            // Position of the first chromosome where the fusion occurs
            BasePair pos_1;
            
            // Position of the second chromosome where the fusion occurs
            BasePair pos_2;
            
            // Orientation for the first chromosome
            Orientation ori_1;
            
            // Orientation for the second chromosome
            Orientation ori_2;
            
            // Number of reads that span the fusiom
            Coverage dd;
            
        };

        typedef std::function<void(const Fusion &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif