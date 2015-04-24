#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "locus.hpp"
#include "feature.hpp"

namespace Spike
{
    enum Mixture
    {
        MixA,
        MixB
    };

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
        
        // Fold ratio from r to v
        Fold fold;

        // Each mixture represents a transcript for a gene
        GeneID geneID;

        // Reference and variant mixtures
        Sequin r, v;
    };

    typedef std::map<SequinID, Sequins> SequinsMap;
    typedef std::map<TranscriptID, Sequin> SequinMap;
    
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
             * RNA data
             */

            SequinsMap r_seqs_gA;
            SequinsMap r_seqs_gB;

            SequinMap r_seqs_iA;
            SequinMap r_seqs_iB;

            // Reference genes
            std::vector<Gene> genes;

            // Reference features
            std::vector<Feature> fs;

            // Reference exons
            std::vector<Feature> exons;

            // Reference introns (spliced junctions)
            std::vector<Feature> introns;

            std::map<TranscriptID, GeneID> iso2Gene;
        
            /*
             * DNA data
             */

            std::map<GeneID, Variation> d_vars;

            /*
             * Metagenomic data
             */

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