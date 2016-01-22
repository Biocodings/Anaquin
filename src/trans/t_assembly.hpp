#ifndef T_ASSEMBLY_HPP
#define T_ASSEMBLY_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TAssembly : Analyzer
    {
        struct Options : FuzzyOptions
        {
            // Reference for the sequins
            FileName chrT;

            // Reference for the endogenous
            FileName endo;
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
                 * Intron-chain statistics
                 */
                
                double cSP, cSN, cFSP, cFSN;
                
                /*
                 * Transcript statistics
                 */
                
                double tSP, tSN, tFSP, tFSN;
                
                /*
                 * Missing statistics
                 */
                
                double mExonN,   mExonR,   mExonP;
                double mIntronN, mIntronR, mIntronP;
                
                /*
                 * Novel statistics
                 */
                
                double nExonN,   nExonR,   nExonP;
                double nIntronN, nIntronR, nIntronP;
            };
            
            std::map<ChromoID, Data> data;
            
            // Reference file for each chromosome
            std::map<ChromoID, FileName> refs;
            
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

        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif