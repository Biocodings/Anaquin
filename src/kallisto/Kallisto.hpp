#ifndef KALLISTO_HPP

#include <map>
#include <set>

namespace Anaquin
{
    struct KMVariant
    {
        // Reference k-mers for reference standard
        std::set<std::string> R;
        
        // Reference k-mers for variant standard
        std::set<std::string> V;
    };
    
    struct KMStats
    {
        // Number of reads estimated to be sequins
        unsigned nSeq = 0;
        
        // Number of reads estimated to be genome (not sequins)
        unsigned nGen = 0;
        
        // Eg: List of reference k-mers spanning variants
        std::map<std::string, KMVariant> vars;
        
        std::map<std::string, unsigned> spans;
        
#ifdef DEBUG
        // All k-mers (for debugging)
        std::map<std::string, unsigned> all;
#endif
    };
}

#endif
