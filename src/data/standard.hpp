#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "data/region.hpp"
#include "data/reader.hpp"
#include "data/sequin.hpp"
#include "data/feature.hpp"
#include "parsers/parser_bed.hpp"

namespace Anaquin
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
            typedef std::map<BaseID, Base_> BaseMap;
            typedef std::map<SequinID, Sequin> SequinMap;

            static Standard& instance(bool reload = false)
            {
                static Standard s;
                
                // Reload the default resources
                if (reload)
                {
                    s = Standard();
                }
                
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

            // The name of the chromosome
            ChromoID id = "chrT";

            // The location of the chromosome
            Locus l;

            /*
             * Shared variables
             */

            // Set of prefix IDs shared between all tools
            std::set<BaseID> baseIDs;
        
            // Set of sequinIDs shared between all tools
            std::set<SequinID> seqIDs;
        
            // Mapping between sequins to the base shared between all tools
            std::map<SequinID, BaseID> seq2base;

            // Sequins for first and second sample
            SequinMap seqs_1, seqs_2;

            // Bases for first and second sample
            BaseMap bases_1, bases_2;

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

            void r_ref(const Reader &);
            void r_mix(const Reader &);

            /*
             * Variant data
             */

            inline const SequinMap &v_seq(Mixture mix) const
            {
                return mix == MixA ? seqs_1 : seqs_1;
            }

            void v_ref(const Reader &);
            void v_mix(const Reader &);

            typedef std::map<BaseID, VariantBase> VariantBaseMap;

            // Bases for mixture A and B
            VariantBaseMap v_seqs_bA, v_seqs_bB;

            // Indexed by locus
            std::map<Locus, Variation> v_vars;

            /*
             * Ladder data
             */

            void l_mix(const Reader &);

            /*
             * Fusion data
             */

            void f_ref(const Reader &);
            void f_mix(const Reader &);

            // Known fusion break-points
            std::set<FusionBreak> f_breaks;

            /*
             * Metagenomic data
             */

            void m_ref(const Reader &);
            void m_mix(const Reader &);

            // Sequin IDs for metagenomic standards
            std::set<SequinID> m_seqIDs;

            // Metagenomic annotation
            std::vector<BedFeature> m_model;

        private:
            Standard();
            Standard(Standard const&) = delete;
    };
}

#endif