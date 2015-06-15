#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "data/reader.hpp"
#include "data/sequin.hpp"
#include "data/feature.hpp"
#include "parsers/parser_bed.hpp"

namespace Spike
{
    #define CHECK_AND_SORT(t) { assert(!t.empty()); std::sort(t.begin(), t.end(), [](const Feature& x, const Feature& y) { return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end); }); }

    enum Mixture
    {
        MixA,
        MixB
    };

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

            // The name of the chromosome (should be chrT)
            ChromoID id;

            // The location of the chromosome
            Locus l;

            /*
             * RNA data
             */

            SequinMap r_seqs_iA, r_seqs_iB;
            PairMap r_seqs_gA, r_seqs_gB;

            std::vector<Sequin>  r_sequins;
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

            // DNA sequins
            SequinMap d_seq_A, d_seq_B;

            // DNA pairs
            PairMap d_pair_A, d_pair_B;

            /*
             * Metagenomic data
             */

            void meta_mod(const Reader &);
            void meta_mix(const Reader &);

            // Metagenomic sequins
            SequinMap m_seq_A, m_seq_B;

            // Metagenomic annotation
            std::vector<BedFeature> m_model;

        private:
            Standard();
            Standard(Standard const&)       = delete;
            void operator=(Standard const&) = delete;

            // Apply default resources for RNA
            void rna();

            // Apply default resources for DNA
            void dna();

            // Apply default resources for metagenomics
            void meta();
    };
}

#endif