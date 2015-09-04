#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "data/reader.hpp"
#include "data/sequin.hpp"
#include "data/feature.hpp"
#include "data/reference.hpp"
#include "parsers/parser_bed.hpp"

namespace Anaquin
{
    #define CHECK_AND_SORT(t) { assert(!t.empty()); std::sort(t.begin(), t.end(), [](const Feature& x, const Feature& y) { return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end); }); }

    struct Variation
    {
        operator const Locus &() const { return l; }
        
        inline bool operator<(const Locus &x) const { return l < x; }

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
            typedef std::vector<Feature>       Features;
            typedef std::map<BaseID, BaseSeq>  BaseMap;
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

            // The name of the chromosome
            ChromoID id = "chrT";

            // The location of the chromosome
            //Locus l;

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

            // Primary features
            Features fs_1;

            /*
             * RNA data
             */

            void r_ref(const Reader &);
            void r_mix(const Reader &);

            TransReference r_trans;

            /*
             * Variant data
             */

            void v_mix(const Reader &);
        
            // Apply reference for known variants
            void v_var(const Reader &);

            // Apply reference for variant standard
            void v_std(const Reader &);

            // Indexed by locus
            std::map<Locus, Variation> v_vars;

            std::set<Variation> __v_vars__;

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

            void m_mix_1(const Reader &);
            void m_mix_2(const Reader &);

            Reference<> r_meta;

        private:
            Standard();
            Standard(Standard const&) = delete;
    };
}

#endif