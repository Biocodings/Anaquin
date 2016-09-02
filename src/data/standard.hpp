#ifndef STANDARD_HPP
#define STANDARD_HPP

#include <set>
#include <memory>
#include "data/reader.hpp"
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

            RnaRef r_trans;

            /*
             * ---------------- Variant analysis ----------------
             */

            // Add reference mixture for VarQuin
            void addVMix(const Reader &);
        
            // Add reference variants for VarQuin (synthetic and genomic/user)
            void addVVar(const Reader &);

            // Add reference standards for VarQuin
            void addVStd(const Reader &);

            VarRef r_var;

        private:
            Standard() {}
            Standard(Standard const&) = delete;
    };
}

#endif