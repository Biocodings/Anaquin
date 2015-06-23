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
        SequinID id;
        
        // The reference position, with the 1st base having position 1
        Locus l;

        // Type of the mutation
        Mutation type;

        Sequence ref, alt;

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
            typedef std::map<BaseID, Base> BaseMap;
            typedef std::map<SequinID, Sequin> SequinMap;

            static Standard& instance()
            {
                static Standard s;
                return s;
            }

            inline const SequinMap& r_sequin(Mixture mix) const
            {
                return mix == MixA ? r_seqs_A : r_seqs_B;
            }

            inline const BaseMap & r_gene(Mixture mix) const
            {
                return mix == Mixture::MixA ? r_seqs_gA : r_seqs_gB;
            }

            // The name of the chromosome (should be chrT)
            ChromoID id;

            // The location of the chromosome
            Locus l;

            /*
             * RNA data
             */

            // Sequins for mixture A and B
            SequinMap r_seqs_A,  r_seqs_B;

            // Genes for mixture A and B
            BaseMap r_seqs_gA, r_seqs_gB;

            // Unique set of sequin names
            std::set<SequinID> r_sequinIDs;
        
            // Sequins and their positions
            std::map<SequinID, Locus> r_sequins;

            std::vector<Feature> r_genes;
            std::vector<Feature> r_exons;
            std::vector<Feature> r_introns;

            std::vector<RNALocus> r_l_exons;
            BasePair r_c_exons;

            // Mapping from sequin to gene
            std::map<SequinID, GeneID> r_isoformToGene;

            void rna_mod(const Reader &);
            void rna_mix(const Reader &);

            /*
             * DNA data
             */

            inline const SequinMap &d_seq(Mixture mix) const
            {
                return mix == MixA ? d_seqs_A : d_seqs_B;
            }

            void dna_mod(const Reader &);
            void dna_mix(const Reader &);

            // Sequins for mixture A and B
            SequinMap d_seqs_A, d_seqs_B;
        
            // Pairs for mixture A and B
            BaseMap d_seqs_bA, d_seqs_bB;

            // Indexed by the position
            std::map<Locus, Variation> d_vars;

            // Unique set of sequin names
            std::set<SequinID> d_sequinIDs;

            /*
             * Metagenomic data
             */

            void meta_mod(const Reader &);
            void meta_mix(const Reader &);

            // Sequins for mixture A and B
            SequinMap m_seqs_A, m_seqs_B;

            // Pairs for mixture A and B
            BaseMap m_seqs_bA, m_seqs_bB;

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