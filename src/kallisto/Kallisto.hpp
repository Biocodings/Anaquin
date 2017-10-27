#ifndef KALLISTO_HPP

namespace Anaquin
{
    struct KMStats
    {
        // Number of reads estimated to be sequins
        unsigned nSeq = 0;
        
        // Number of reads estimated to be genome (not sequins)
        unsigned nGen = 0;
        
#ifdef DEBUG
        // All k-mers (for debugging)
        std::map<std::string, unsigned> all;
#endif
    };
    
    KMStats Kallisto(const std::string &, const std::string &, const std::string &);
}

#endif
