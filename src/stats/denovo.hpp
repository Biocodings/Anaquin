#ifndef GI_DENOVO_HPP
#define GI_DENOVO_HPP

#include <vector>
#include "types.hpp"

namespace Spike
{
    struct Contig
    {
        BasePair l;
        std::string id;
        std::string seq;
    };

    struct DNStats
    {
        std::vector<Contig> contigs;

        BasePair min, max;
        BasePair mean, sum;
        BasePair N20, N50, N80;
    };

    struct DNAsssembly
    {
        static DNStats stats(const std::string &file);
    };
}

#endif