#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "locus.hpp"
#include "feature.hpp"

namespace Spike
{
    enum Group { A, B, C, D };

    enum Zygosity
    {
        Homozygous,
        Heterzygous,
    };
    
    typedef std::string Sequence;

    struct Variation
    {
        GeneID id;        
        BasePair pos;
        Zygosity zy;
        Sequence r, m;
    };
    
    struct Gene
    {
        GeneID id;
        
        // 1-indexed locus relative to the chromosome
        Locus l;

        std::vector<Feature> exons;
        std::vector<Feature> introns;
    };
    
    struct Sequin
    {
        IsoformID id;        
        Locus l;

        // Fold ratio relative to the other sequin
        Fold fold;

        // Amount of abundance
        Concentration raw;

        // Amount of abundance after normalization
        Concentration fpkm;
    };

    struct Sequins
    {
        inline Concentration raw() const  { return r.raw + v.raw;   }
        inline Concentration fpkm() const { return r.fpkm + v.fpkm; }

        Group grp;

        // Each mixture represents a transcript for a gene
        GeneID geneID;

        // Reference and variant mixtures
        Sequin r, v;
    };

    class Standard
    {
        public:
            static Standard& instance()
            {
                static Standard s;
                return s;
            }

            ChromoID id;
        
            // The location of the chromosome
            Locus l;
        
            /*
             * Structural variations
             */

            std::vector<Variation> vars;

            /*
             * Metagenomic sequins
             */

            /*
             * RNA sequins
             */

            std::map<GeneID, Sequins> seqs_gA;
            std::map<GeneID, Sequins> seqs_gB;

            std::map<TranscriptID, Sequin> seqs_iA;
            std::map<TranscriptID, Sequin> seqs_iB;
        
            std::map<TranscriptID, GeneID> iso2Gene;
        
            // Reference genes
            std::vector<Gene> genes;
        
            // Reference features
            std::vector<Feature> fs;

            // Reference exons
            std::vector<Feature> exons;
        
            // Reference introns (spliced junctions)
            std::vector<Feature> introns;

        private:
            Standard();
            Standard(Standard const&)       = delete;
            void operator=(Standard const&) = delete;
    };
}

#endif