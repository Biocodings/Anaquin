#ifndef STANDARD_HPP
#define STANDARD_HPP

#include <memory>
#include "data/reader.hpp"
#include "data/feature.hpp"
#include "data/reference.hpp"

namespace Anaquin
{
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

            static constexpr const char * chrT = "chrT";

            /*
             * ---------------- Transcriptome analysis ----------------
             */

            // Add a reference annotation
            void addTRef(const Reader &);

            void addTMix(const Reader &);
            void addTDMix(const Reader &);

            TransRef r_trans;

            /*
             * ---------------- Variant analysis ----------------
             */

            void addVMix(const Reader &);
        
            // Add known variants
            void addVar(const Reader &);

            // Add standards to the reference
            void addStd(const Reader &);

            // Add intervals to the reference
            void addInters(const Reader &);
        
            VarRef r_var;
        
            /*
             * ---------------- Ladder analysis ----------------
             */

            void addLMix(const Reader &);

            LadderRef r_lad;

            /*
             * ---------------- Fusion analysis ----------------
             */
        
            // Fusions for FuseQuin
            void addFRef(const Reader &);

            // Standards for FuseQuin
            void addFStd(const Reader &);

            // Splicing for FuseQuin
            void addFSplice(const Reader &);
        
            // Mixture for FuseQuin
            void addFMix(const Reader &);

            FusionRef r_fus;
        
            /*
             * ---------------- Metagenomic analysis ----------------
             */

            void addMRef(const Reader &);
            void addMMix(const Reader &);

            MetaRef r_meta;

        private:
            Standard() {}
            Standard(Standard const&) = delete;
    };
}

#endif