#ifndef V_ALIGN_HPP
#define V_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        typedef AnalyzerOptions Options;
        
        struct Stats : public AlignmentStats
        {
            // Returns the overall sensitivity
            inline Proportion sn(const ChromoID &cID) const
            {
                return static_cast<Proportion>(covered(cID)) / length(cID);
            }

            // Returns the individual sensitivity
            inline Proportion sn(const ChromoID &cID, const SequinID &sID) const
            {
                // How many bases covered?
                const auto covered = data.at(cID).covered.at(sID);
                
                // How long is the reference?
                const auto length  = data.at(cID).length.at(sID);

                return static_cast<Proportion>(covered) / length;
            }

            // Returns the overall precision
            inline Proportion pc(const ChromoID &cID) const
            {
                const auto &tp = data.at(cID).tp;
                const auto &fp = data.at(cID).fp;

                return static_cast<Proportion>(tp) / (tp + fp);
            }

            inline Base length(const ChromoID &cID)  const { return sum(data.at(cID).length);  }
            inline Base covered(const ChromoID &cID) const { return sum(data.at(cID).covered); }
            
            struct Data
            {
                Counts tp, fp;
                
                // Number of bases covered
                std::map<SequinID, Base> covered;
                
                // Number of bases for each reference
                std::map<SequinID, Base> length;
                
                // Distribution for the sequins
                SequinHist hist;                
            };
            
            std::map<ChromoID, Data> data;

            // Absolute detection limit
            Limit limit;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif