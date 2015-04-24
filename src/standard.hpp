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
        inline Concentration abund(bool norm) const
        {
            return norm ? fpkm : raw;
        }

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
        inline Concentration abund(bool norm) const
        {
            return r.abund(norm) + v.abund(norm);
        }
        
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
             * RNA sequins
             */

            std::map<GeneID, Sequins> r_seqs_gA;
            std::map<GeneID, Sequins> r_seqs_gB;

            std::map<TranscriptID, Sequin> r_seqs_iA;
            std::map<TranscriptID, Sequin> r_seqs_iB;
        
            /*
             * DNA sequins
             */
        
            std::vector<Variation> vars;
        
            /*
             * Metagenomic sequins
             */

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

            void rna();
            void dna();
            void meta();
    };
}

#endif