#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "data/region.hpp"
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
        operator const Locus &() const { return l; }

        SequinID id;
        
        // The reference position, with the 1st base having position 1
        Locus l;

        // Type of the mutation
        Mutation type;

        Sequence ref, alt;

        Genotype gt;
        
        // Allelle frequency
        double af;

        // Allele count in genotypes
        Counts ac;
        
        // Total number of alleles in called genotypes
        Counts an;
        
        // Combined depth across samples
        unsigned dp;
        
        // Depth for reference
        unsigned dp_r;
        
        // Depth for alternative
        unsigned dp_a;
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

            // Sequin IDs for RNA standards
            std::set<SequinID> r_seqIDs;

            // Sequins and their positions
            std::map<SequinID, Locus> r_sequins;

            std::vector<Feature> r_exons;
            std::vector<Feature> r_genes;
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

            typedef std::map<BaseID, VariantBase> VariantBaseMap;

            // Bases for mixture A and B
            VariantBaseMap d_seqs_bA, d_seqs_bB;

            // Indexed by the position
            std::map<Locus, Variation> d_vars;

            // Sequin IDs for DNA standards
            std::set<SequinID> d_seqIDs;

            /*
             * Ladder data
             */

            void ladder_mix(const Reader &);

            // Sequins for mixture A and B
            SequinMap l_seqs_A, l_seqs_B;

            std::map<TypeID, SequinID> l_map;

            /*
             * Fusion data
             */

            void f_mod(const Reader &);
            void f_mix(const Reader &);

            // Locus for each fusion in the forward strand
            std::map<SequinID, Locus> f_f_fusions;

            // Locus for each fusion in the reverse strand
            std::map<SequinID, Locus> f_r_fusions;

            // Mixture for each fusion sequin
            SequinMap f_seqs_A;

            /*
             * Metagenomic data
             */

            void meta_mod(const Reader &);
            void meta_mix(const Reader &);

            // Sequins for mixture A and B
            SequinMap m_seqs_A, m_seqs_B;

            // Bases for mixture A and B
            BaseMap m_seqs_bA, m_seqs_bB;

            // Sequin IDs for metagenomic standards
            std::set<SequinID> m_seqIDs;

            // Metagenomic annotation
            std::vector<BedFeature> m_model;

        private:
            Standard();
            Standard(Standard const&)       = delete;
            void operator=(Standard const&) = delete;

            // Apply resources for RNA
            void rna();

            // Apply resources for variant
            void variant();

            // Apply resources for metagenomics
            void meta();

            // Apply resources for conjoint
            void ladder();

            // Apply resources for fusion
            void fusion();

            // Apply resources for clinical
            void clinical();

            // Apply resources for cancer
            void cancer();
    };
}

#endif