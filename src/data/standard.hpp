#ifndef GI_STANDARD_HPP
#define GI_STANDARD_HPP

#include <map>
#include <vector>
#include <memory>
#include "data/reader.hpp"
#include "data/feature.hpp"
#include "data/reference.hpp"
#include "parsers/parser_bed.hpp"

namespace Anaquin
{
    #define CHECK_AND_SORT(t) { assert(!t.empty()); std::sort(t.begin(), t.end(), [](const Feature& x, const Feature& y) { return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end); }); }

    class Standard
    {
        public:
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

            VarRef r_var;
        
            /*
             * Ladder data
             */

            void l_mix(const Reader &);

            LadderRef r_lad;
        
            /*
             * Fusion data
             */

            void f_ref(const Reader &);
            void f_mix(const Reader &);

            FusionRef r_fus;
        
            /*
             * Metagenomic data
             */

            void m_ref(const Reader &);
            void m_mix(const Reader &);

            MetaRef r_meta;

        private:
            Standard() {}
            Standard(Standard const&) = delete;
    };
}

#endif