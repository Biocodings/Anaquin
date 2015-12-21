#ifndef T_ASSEMBLY_HPP
#define T_ASSEMBLY_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TAssembly : Analyzer
    {
        struct Options : FuzzyOptions
        {
            // Path for the reference and query GTF
            FileName ref, query;
        };

        struct Stats : public MappingStats
        {
            struct Data
            {
                /*
                 * Base statistics
                 */
                
                double bSP, bSN;
                
                /*
                 * Exon statistics
                 */
                
                double eSP, eSN, eFSP, eFSN;

                /*
                 * Intron statistics
                 */
                
                double iSP, iSN, iFSP, iFSN;

                /*
                 * Transcript statistics
                 */
                
                double tSP, tSN, tFSP, tFSN;
            };
            
            std::map<ChromoID, Data> data;
            
            /*
             * Statistics for detection limit (chrT only)
             */
            
            Limit bLimit, eLimit, iLimit, tLimit;
            SequinHist bHist, eHist, iHist, tHist;
        };

        // Analyze a single sample
        static Stats analyze(const FileName &, const Options &o = Options());
        
        // Analyze multiple replicates
        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o = Options())
        {
            std::vector<TAssembly::Stats> stats;
            
            for (auto &file : files)
            {
                stats.push_back(analyze(file, o));
            }
            
            return stats;
        }

        static void report(const FileName &, const Options &o = Options());
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif