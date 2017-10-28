#ifndef KALLISTO_HPP

#include <map>
#include <vector>
#include "data/data.hpp"

namespace Anaquin
{
    struct KMPair
    {
        Kmer normal, revComp;
    };
    
    struct KMVariant
    {
        // Reference k-mers for reference standard
        std::vector<KMPair> R;
        
        // Reference k-mers for variant standard
        std::vector<KMPair> V;
    };
    
    struct KMStats
    {
        // Number of reads estimated to be sequins
        unsigned nSeq = 0;
        
        // Number of reads estimated to be genome (not sequins)
        unsigned nGen = 0;
        
        // Eg: List of reference k-mers spanning variants
        std::map<std::string, KMVariant> vars;
        
        // Measured counts for reference spanning k-mers
        std::map<std::string, unsigned> spans;
        
#ifdef DEBUG
        // All k-mers (for debugging)
        std::map<std::string, unsigned> all;
#endif
    };
}

#endif
