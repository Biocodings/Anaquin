#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "locus.hpp"
#include "feature.hpp"

namespace Spike
{
    #define CHECK_AND_SORT(t) { assert(!t.empty()); std::sort(t.begin(), t.end(), [](const Feature& x, const Feature& y) { return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end); }); }

    enum Mixture
    {
        MixA,
        MixB
    };

    enum Group { A, B, C, D };

    typedef std::string Sequence;

    struct Variation
    {
        ChromoID id;
        
        // The reference position, with the 1st base having position 1
        Locus l;
        
        // Type of the mutation
        Mutation m;
        
        Sequence r, a;
        
        Genotype gt;
        
        // Allelle frequency
        Counts af;
        
        // Allele count in genotypes
        Counts ac;
        
        // Total number of alleles in called genotypes
        Counts an;
        
        // Combined depth across samples
        unsigned dp;
    };
    
    struct Sequin
    {
        operator Locus()     const { return l;  }
        operator IsoformID() const { return id; }

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
            typedef std::map<SequinID, Sequin>  SequinMap;
            typedef std::map<SequinID, Sequins> SequinsMap;

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

            std::vector<Sequin> r_sequins;

            std::vector<Feature> r_fs;
            std::vector<Feature> r_genes;
            std::vector<Feature> r_exons;
            std::vector<Feature> r_introns;

            std::vector<RNALocus> r_l_exons;
            BasePair r_c_exons;

            std::map<TranscriptID, GeneID> r_iso2Gene;

            /*
             * DNA data
             */

            std::set<std::string> d_seqs;
            std::map<Locus, Variation> d_vars;

            /*
             * Metagenomic data
             */

        private:
            Standard();
            Standard(Standard const&)       = delete;
            void operator=(Standard const&) = delete;

            void rna (const std::string &mix = "data/rna/rna_mixtures.csv");
            void dna (const std::string &mix = "data/rna/rna_mixtures.csv");
            void meta(const std::string &mix = "data/rna/meta_mixtures.csv");
    };
}

#endif