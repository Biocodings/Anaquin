#ifndef V_KSTATS_HPP
#define V_KSTATS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VKStats
    {
        struct Stats
        {
            // Number of reads estimated to be sequins
            unsigned nSeq = 0;
            
            // Number of reads estimated to be genome (not sequins)
            unsigned nGen = 0;

            inline float dilution() const
            {
                return (float) nSeq / (nSeq + nGen);
            }
            
            // Ladder for allelle frequency
            SequinStats af;
        };
        
        struct Options : public AnalyzerOptions
        {
            // Reference FASTA file for sequins
            FileName sFA;
        };
        
        static Stats analyze(const std::vector<FileName> &, const Options &o);
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
