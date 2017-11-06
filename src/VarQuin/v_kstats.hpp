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
                // Number of reads mapping for each sequin
                std::map<SequinID, Counts> s2r;
                
                // Estimated sequin abundance
                std::map<SequinID, Coverage> s2a;
                
                // Raw counts
                std::map<SequinID, std::vector<Counts>> raws;
                
                // Minimum, maximum, medians and standard deviation
                std::map<SequinID, Counts> mins, maxs, meds, sds;
            };
            
            Abundance R, F;
            
            // Length for sequins
            std::map<SequinID, Base> s2l;
            
            // Counts for ssequin
            std::map<SequinID, std::vector<Counts>> s2c;
            
            inline float dilution() const
            {
                return (float) kStats.R.nMatch / (kStats.R.nMatch + kStats.R.nNMatch);
            }
            
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
