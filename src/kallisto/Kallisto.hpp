#ifndef KALLISTO_HPP
#define KALLISTO_HPP

#include <set>
#include <map>
#include <vector>
#include "data/data.hpp"

namespace Anaquin
{
    typedef std::set<SequinID> Matches;
    typedef std::map<SequinID, std::map<Kmer, Counts>> SequinCounts;
    
    struct KStats
    {
        // Sequin standards (no _R and _V)
        std::set<SequinID> seqs;
        
        struct KAbund
        {
            // Measured counts for unique k-mers
            SequinCounts uniqs;
            
            // Measured counts for shared k-mers
            SequinCounts shared;

            // Counts for sequins
            std::map<SequinID, Counts> m;

            // Number of reads matching this index (e.g. sequins)
            Counts nMatch = 0;
            
            // Number of reads not matching this index (e.g. genome)
            Counts nNMatch = 0;
        };
        
        // Abundance for the reverse and forward genome
        KAbund R, F;
        
        // Index for all reference sequin k-mers
        FileName i1;
    };

    void KInit(const FileName &, const FileName &, unsigned);
    
    FileName KHumanFA(const FileName &, std::map<SequinID, Base> &);

    // Build Kallisto index from a FASTA file
    FileName KBuildIndex(const FileName &, unsigned k = 31);
    
    // Heavily modified Kallisto k-mer counting
    KStats KCount(const FileName &, const FileName &, const FileName &, const FileName &, unsigned k = 31);

    // Query a k-mer from a Kallisto index
    Matches KQuery(const FileName &, const Sequence &s, unsigned k = 31);
}

#endif
