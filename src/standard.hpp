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

    struct TGene
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

    class Standard
    {
        public:
            typedef std::map<SequinID, Sequins> SequinsMap;
            typedef std::map<TranscriptID, Sequin> SequinMap;

            static Standard& instance()
            {
                static Standard s;
                return s;
            }

            inline const SequinsMap& r_pair(Mixture mix) const
            {
                return mix == MixA ? r_seqs_gA : r_seqs_gB;
            }

            inline const SequinMap& r_sequin(Mixture mix) const
            {
                return mix == MixA ? r_seqs_iA : r_seqs_iB;
            }

            ChromoID id;

            // The location of the chromosome
            Locus l;

            /*
             * RNA data
             */

            SequinMap  r_seqs_iA, r_seqs_iB;
            SequinsMap r_seqs_gA, r_seqs_gB;
        
            std::vector<Feature> r_fs;
            std::vector<TGene>   r_genes;
            std::vector<Feature> r_exons;
            std::vector<Feature> r_introns;

            std::map<TranscriptID, GeneID> r_iso2Gene;

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

            void rna(const std::string &mix = "data/rna/rna_mixtures.csv");
            void dna();
            void meta();
    };
}

#endif