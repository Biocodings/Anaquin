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

            // Reference genomic chromosomes
            static std::set<ChrID> genoIDs;
        
            /*
             * ---------------- General analysis ----------------
             */

            // Whether the chromosome is synthetic
            static bool isSynthetic(const ChrID &);

            // Whether the chromosome is genomic
            static bool isGenomic(const ChrID &);
        
            // Add genomic chromosome, which can be checked by isGenomic()
            static void addGenomic(const ChrID &);
        
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

            // Add reference mixture for VarQuin
            void addVMix(const Reader &);
        
            // Add reference variants for VarQuin (synthetic and genomic/user)
            void addVVar(const Reader &);

            // Add reference standards for VarQuin
            void addVStd(const Reader &);

            // Add intervals to the reference
            void addInters(const Reader &);

            VarRef r_var;

            /*
             * ---------------- Ladder analysis ----------------
             */

            void addLMix(const Reader &);

            LadderRef r_lad;

            /*
             * ---------------- Structual analysis ----------------
             */
        
            void addSStruct(const Reader &);

            StructRef r_str;
        
            /*
             * ---------------- Fusion analysis ----------------
             */
        
            // Fusions for FuseQuin
            void addFRef(const Reader &);

            // Standards for FuseQuin
            void addFStd(const Reader &);

            // Intron junctions for FuseQuin
            void addFJunct(const Reader &);
        
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