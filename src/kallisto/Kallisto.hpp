#ifndef KALLISTO_HPP
#define KALLISTO_HPP

#include <map>
#include <vector>
#include "data/data.hpp"

namespace Anaquin
{
    struct KMPair
    {
        Kmer norm, rcom;
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
        // Index for all reference sequin k-mers
        FileName i1;
        
        // Number of reads estimated to be sequins
        unsigned nSeq = 0;
        
        // Number of reads estimated to be genome (not sequins)
        unsigned nGen = 0;
        
        // Eg: List of reference k-mers spanning variants
        std::map<SequinID, KMVariant> vars;
        
        // Measured counts for reference spanning k-mers
        std::map<Kmer, unsigned> spans;
        
        // Measured counts for all reference sequin k-mers
        std::map<Kmer, Counts> k2c;
    };

    // Query a k-mer for all reference sequin k-mers
    bool KQuerySeqs(const Sequence &, unsigned);

    // Query a k-mer from a Kallisto index
    bool KQuery___(const FileName &, const Sequence &);

    FileName KHumanFASTA(const FileName &file);
    
    // Translate k-mer to sequin
    SequinID KKM2Sequin(const Kmer &, unsigned);
    
    // Build Kallisto index from a FASTA file
    FileName KBuildIndex(const FileName &, unsigned k);

    // Heavily modified Kallisto k-mer counting
    KMStats KCount(const FileName &, const FileName &, const FileName &, const FileName &, unsigned);
}

#endif
