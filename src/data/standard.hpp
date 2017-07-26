#ifndef STANDARD_HPP
#define STANDARD_HPP

#include "data/ladder.hpp"
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

            // Add sequin regions in BED format
            static BedData readBED(const Reader &, Base trim = 0);

            // Add sequin regions in GTF format
            static std::shared_ptr<GTFData> readGTF(const Reader &);

            Ladder readLength(const Reader &);
        
           /*
             * ---------------- RnaQuin analysis ----------------
             */

            Ladder readGene(const Reader &);
            Ladder readIsoform(const Reader &);
            Ladder readGeneL(const Reader &);
        
            RnaRef r_rna;

            /*
             * ---------------- MetaQuin analysis ----------------
             */

            Ladder addMMix(const Reader &);
        
            MetaRef r_meta;
        
            /*
             * ---------------- VarQuin analysis ----------------
             */

            // Add reference variants in VCF format
            inline void addVVar(const Reader &r) { r_var.readVRef(r); }

            Ladder addAF(const Reader &);
            Ladder addCNV(const Reader &);
            Ladder addCon1(const Reader &);
            Ladder addCon2(const Reader &);
        
            VarRef r_var;

        private:
            Standard() {}
            Standard(Standard const&) = delete;
    };
}

#endif
