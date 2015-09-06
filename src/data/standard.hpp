#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include "data/reader.hpp"
#include "data/sequin.hpp"
#include "data/feature.hpp"
#include "data/reference.hpp"
#include "data/variation.hpp"
#include "parsers/parser_bed.hpp"

namespace Anaquin
{
    #define CHECK_AND_SORT(t) { assert(!t.empty()); std::sort(t.begin(), t.end(), [](const Feature& x, const Feature& y) { return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end); }); }

    class Standard
    {
        public:
            typedef std::vector<Feature>       Features;
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

            /*
             * RNA data
             */

            void r_ref(const Reader &);
            void r_mix(const Reader &);

            TransRef r_trans;

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

            VarRef v_ref;
        
            /*
             * Ladder data
             */

            void l_mix(const Reader &);

            LadderRef l_ref;
        
            /*
             * Fusion data
             */

            void f_ref(const Reader &);
            void f_mix(const Reader &);

            // Known fusion break-points
            std::set<FusionBreak> f_breaks;

            FusionRef r_fus;
        
            /*
             * Metagenomic data
             */

            void m_mix_1(const Reader &);
            void m_mix_2(const Reader &);

            MetaRef r_meta;

        private:
            Standard() {}
            Standard(Standard const&) = delete;
    };
}

#endif