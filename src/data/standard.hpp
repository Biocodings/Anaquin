#ifndef STANDARD_HPP
#define STANDARD_HPP

#include "data/vData.hpp"
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

            Ladder readGDiff(const Reader &);
            Ladder readIDiff(const Reader &);

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

            static VCFLadder addVCF(const Reader &r);

            Ladder addAF(const Reader &);   // From mixture file
            Ladder addCNV(const Reader &);  // From mixture file
            Ladder addCon1(const Reader &); // From mixture file
            Ladder addCon2(const Reader &); // From mixture file
        
            Translate addSeq2Unit(const Reader &);
            Translate addUnit2Seq(const Reader &);

            VarRef r_var;

        private:
            Standard() {}
            Standard(Standard const&) = delete;
    };
}

#endif
