#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "locus.hpp"
#include "feature.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
    enum Group { A, B, C, D };
    
    struct Variation
    {
        GeneID id;
        Locus l;
        std::string r;
        std::string m;
    };
    
    struct Gene
    {
        GeneID id;
        
        // Location of the gene relative to the chromosome
        Locus l;
        
        std::vector<Feature> exons;
        std::vector<Feature> introns;
    };
    
    struct Sequin
    {
        TranscriptID id;

        // Fold ratio relative to the other mixture
        Fold fold;

        // Amount of abundance for this mixture
        Concentration reads;
    };

    struct Sequins
    {
        Group grp;

        // Each mixture represents a transcript for a gene
        GeneID geneID;

        // Reference and variant mixtures
        Sequin r, v;
    };
    
    struct Standard
    {
        ChromoID id;
        
        // The location of the chromosome
        Locus l;

        std::vector<Variation> vars;
        
        std::map<GeneID, Sequins> seqs_gA;
        std::map<GeneID, Sequins> seqs_gB;

        std::map<TranscriptID, Sequin> seqs_iA;
        std::map<TranscriptID, Sequin> seqs_iB;

        std::map<TranscriptID, GeneID> iso2Gene;
        
        // Known genes
        std::vector<Gene> genes;
        
        // Known features
        std::vector<Feature> fs;
        
        // Known exons
        std::vector<Feature> exons;
        
        // Known introns (spliced junctions)
        std::vector<Feature> introns;
    };    
}

#endif