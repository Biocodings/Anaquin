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
            // Returns the sensitivity at the base level. This is also the base coverage.
            inline Proportion bSN(const ChromoID &cID, const SequinID &sID) const
            {
                // How many bases covered?
                const auto covered = data.at(cID).covered.at(sID);
                
                // How long is the reference?
                const auto length  = data.at(cID).length.at(sID);

                assert(covered <= length);

                return static_cast<Proportion>(covered) / length;
            }

            // Returns the total number of bases covered for a chromosome
            inline Base covered(const ChromoID &cID) const
            {
                return sum(data.at(cID).covered);
            }

            // Returns the total reference length for a chromosome
            inline Base length(const ChromoID &cID) const
            {
                return sum(data.at(cID).length);
            }
            
            inline Proportion bSN(const ChromoID &cID) const
            {
                return static_cast<Proportion>(covered(cID)) / length(cID);
            }

            struct Data
            {
                Confusion m;
                
                // Number of bases covered
                std::map<SequinID, Base> covered;
                
                // Number of bases for each reference
                std::map<SequinID, Base> length;
            };
            
            std::map<ChromoID, Data> data;

            // Distribution for the sequins
            SequinHist h;

            // Absolute detection limit
            Limit limit;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif