#ifndef V_KSTATS_HPP
#define V_KSTATS_HPP

#include "Kallisto.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VKStats
    {
        struct Stats
        {
            // Raw statistics from Kallisto
            KMStats kStats;
            
            // Counts for each sequin
            std::map<SequinID, std::vector<Counts>> s2c;
            
            inline float dilution() const
            {
                return (float) kStats.nSeq / (kStats.nSeq + kStats.nGen);
            }
            
            // Minimum, maximum and medians
            std::map<SequinID, Counts> mins, maxs, meds;

            // Ladder for allelle frequency
            SequinStats af;
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() : k(31) {}
            
            unsigned k;
            
            // Reference FASTA file for sequins
            FileName fa;
        };
        
        static Stats analyze(const std::vector<FileName> &, const Options &o);
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
