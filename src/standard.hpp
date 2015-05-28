#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "locus.hpp"
#include "feature.hpp"
#include "parsers/parser_bed.hpp"

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

        inline Concentration abund() const { return raw; }

        SequinID id;
        Locus l;
        
        // Amount of abundance
        Concentration raw;
    };

    struct Sequins
    {
        inline Concentration abund() const
        {
            return r.abund() + v.abund();
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
            typedef std::map<SequinID, Sequin>  SequinMap;
            typedef std::map<SequinID, Sequins> PairMap;

            static Standard& instance()
            {
                static Standard s;
                return s;
            }

            inline const PairMap& r_pair(Mixture mix) const
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
            PairMap r_seqs_gA, r_seqs_gB;

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

            inline const PairMap &d_pair(Mixture mix) const
            {
                return mix == MixA ? d_pair_A : d_pair_B;
            }

            inline const SequinMap &d_seq(Mixture mix) const
            {
                return mix == MixA ? d_seq_A : d_seq_B;
            }

            std::map<Locus, Variation> d_vars;

            // DNA annotation
            std::vector<BedFeature> d_annot;

            // DNA sequins for each of the mixture
            SequinMap d_seq_A, d_seq_B;

            // DNA pairs for each of the mixture
            PairMap d_pair_A, d_pair_B;

            /*
             * Metagenomic data
             */

            // Metagenomic sequins for each of the mixture
            SequinMap m_seq_A, m_seq_B;

            // Metagenomic pairs for each of the mixture
            PairMap m_pair_A, m_pair_B;
        
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