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

            // Whether the chromosome is genomic
            static bool isGenomic(const ChrID &cID)
            {
                A_ASSERT(!cID.empty());
                return Standard::genoIDs.count(cID);
            }

            // Add genomic chromosome, which can be checked by isGenomic()
            static void addGenomic(const ChrID &cID)
            {
                A_ASSERT(!cID.empty());
                Standard::genoIDs.insert(cID);
            }

            /*
             * ---------------- RnaQuin analysis ----------------
             */

            // Add reference annotation for RnaQuin in GTF format
            inline void addRRef(const Reader &r) { r_rna.readRef(r); }

            void addRMix(const Reader &);
            void addRDMix(const Reader &);

            RnaRef r_rna;

            /*
             * ---------------- MetaQuin analysis ----------------
             */

            // Add reference mixture for MetaQuin
            void addMMix(const Reader &);

            // Add reference differential mixture for MetaQuin
            void addMDMix(const Reader &);
        
            // Add reference annotation for MetaQuin in BED format
            inline void addMBed(const Reader &r) { r_meta.readBed(r); }

            MetaRef r_meta;
        
            /*
             * ---------------- VarQuin analysis ----------------
             */

            // Add reference variants in VCF format
            inline void addVVar(const Reader &r) { r_var.readVRef(r); }

            // Add sequin regions in BED format
            inline void addVRef(const Reader &r, Base trim = 0) { r_var.readGBRef(r, trim); }
        
            void addCon(const Reader &);
            void addAll(const Reader &);
            void addCNV(const Reader &);
        
            VarRef r_var;

        private:
            Standard() {}
            Standard(Standard const&) = delete;
    };
}

#endif
