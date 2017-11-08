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
            KStats kStats;
            
            struct Abundance
            {
                // Raw counts
                std::map<SequinID, std::vector<Counts>> raws;
                
                // Statistics
                std::map<SequinID, Counts> mins, maxs, meds, sds;
            };
            
            Abundance R, F;
            
            // Length for sequins
            std::map<SequinID, Base> s2l;
            
            inline float dilution() const
            {
                return (float) kStats.R.nMatch / (kStats.R.nMatch + kStats.R.nNMatch);
            }
            
            inline float error() const
            {
                return (float) kStats.R.nNMKMatch / (kStats.R.nNMKMatch + kStats.R.nMKMatch);
            }
            
            // Ladder for allelle frequency
            SequinStats af;
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() : k(31), thr(1) {}
            
            unsigned k;
            
            // Show sequin and genome reads?
            bool showReads;
            
            // Number of threads
            Counts thr;
            
            // Reference FASTA file for sequins
            FileName fa;
        };
        
        static Stats analyze(const std::vector<FileName> &, const Options &o);
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
