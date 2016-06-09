#ifndef V_ALIGN_HPP
#define V_ALIGN_HPP

#include "data/types.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        typedef AnalyzerOptions Options;
        
        struct Stats : public AlignmentStats
        {
            /*
             * Returns the overall sensitivity, defined by the number of bases aligned
             * with respect to the reference.
             */

            inline Proportion sn(const ChrID &cID) const
            {
                return static_cast<Proportion>(covered(cID)) / length(cID);
            }

            // Returns the individual sensitivity
            inline Proportion sn(const ChrID &cID, const SequinID &sID) const
            {
                // How many bases covered?
                const auto covered = data.at(cID).covered.at(sID);
                
                // How long is the reference?
                const auto length  = data.at(cID).length.at(sID);

                return static_cast<Proportion>(covered) / length;
            }

            // Returns the overall precision
            inline Proportion pc(const ChrID &cID) const
            {
                const auto &tp = data.at(cID).tp;
                const auto &fp = data.at(cID).fp;

                return static_cast<Proportion>(tp) / (tp + fp);
            }

            inline Base length(const ChrID &cID)  const { return sum(data.at(cID).length);  }
            inline Base covered(const ChrID &cID) const { return sum(data.at(cID).covered); }
            
            struct Data
            {
                Counts tp, fp;

                // FP alignments
                std::vector<ReadID> afp;

                // Number of bases covered
                std::map<GeneID, Base> covered;

                // Number of bases for each reference
                std::map<GeneID, Base> length;
                
                // Distribution for the genes (eg: sequins)
                Hist hist;
                
                /*
                 * Statistics for genes (eg: sequins)
                 */
                
                std::map<GeneID, Counts> gtp;
                std::map<GeneID, Counts> gfp;
            };

            std::map<ChrID, Data> data;

            /*
             * Genomic statistics
             */
            
            // Genes to covered
            std::map<GeneID, Base> g2c;
            
            // Genes to length
            std::map<GeneID, Base> g2l;

            // Total covered by all genes
            Base gc = 0;
            
            // Total length by all genes
            Base gl = 0;
            
            // Number of TP for genomic genes
            Counts gtp = 0;
            
            // Number of FP for genomic genes
            Counts gfp = 0;
            
            // Overall sensitivity for all genomic genes
            inline Proportion gsn() const
            {
                return static_cast<Proportion>(gc) / (gc + gl);
            }
            
            // Overall precision for all genomic genes
            inline Proportion gpc() const
            {
                return static_cast<Proportion>(gtp) / (gtp + gfp);
            }

            /*
             * Sequin statistics
             */
            
            // Sequins to precision
            std::map<SequinID, Proportion> s2p;
            
            // Sequins to sensitivity
            std::map<SequinID, Proportion> s2s;
            
            // Sequins to covered
            std::map<SequinID, Base> s2c;
            
            // Sequins to length
            std::map<SequinID, Base> s2l;
            
            // Sequins to reads
            std::map<SequinID, Coverage> s2r;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif