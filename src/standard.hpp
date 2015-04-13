#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <set>
#include <map>
#include <list>
#include <vector>
#include "locus.hpp"
#include "feature.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
    enum Group
    {
        A,
        B,
        C,
        D
    };
    
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
    
    struct IMixture
    {
        IsoformID tranID;
        
        // Fold ratio relative to the other mixture
        Fold fold;

        // Amount of abundance for this mixture
        Concentration reads;
    };

    struct GMixture
    {
        Group grp;

        // Each mixture represents a transcript for a gene
        GeneID geneID;

        // Reference and variant mixtures
        IMixture r, v;
    };
    
    struct Standard
    {
        inline bool matchFeature(const Feature &q) const
        {
            for (auto r : fs)
            {
                if (r.type == q.type && r.l.contains(q.l))
                {
                    return true;
                }
            }
            
            return false;
        }
        
        // Whether the given gene is a part of the standard
        bool known(const GeneID &id) const;
        
        ChromoID id;
        
        // The location of the chromosome
        Locus l;
        
        std::vector<Variation> vars;
        
        std::map<GeneID, GMixture> mix_gA;
        std::map<GeneID, GMixture> mix_gB;

        std::map<IsoformID, IMixture> mix_iA;
        std::map<IsoformID, IMixture> mix_iB;

        std::map<IsoformID, GeneID> iso2Gene;
        
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