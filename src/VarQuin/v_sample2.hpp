#ifndef V_SAMPLE2_HPP
#define V_SAMPLE2_HPP

#include "stats/analyzer.hpp"
#include "tools/coverage.hpp"

namespace Anaquin
{
    struct VSample2
    {
        enum class Method
        {
            Mean,
            Median,
            Reads,
            Prop,
        };
        
        struct SampledInfo
        {
            RegionID rID;
            
            // Alignment coverage for the genome
            Coverage gen;
            
            // Alignment coverage before subsampling
            Coverage before;
            
            // Alignment coverage after subsampling
            Coverage after;
            
            // Normalization factor
            Proportion norm;
        };
        
        struct GenomeSequins
        {
            Counts countGen = 0;
            Counts countSyn = 0;
        };
        
        struct Stats
        {
            // Total number of subsampling regions
            Counts count = 0;
            
            // Number of regions without alignment (genomic)
            Counts noGAlign = 0;
            
            // Number of regions without alignment (synthetic)
            Counts noSAlign = 0;
            
            Coverage afterGen,  afterSyn;
            Coverage beforeGen, beforeSyn;
            
            // Summary statistics for normalization
            double normAver, normSD;
            
            GenomeSequins totBefore,  totAfter;
            GenomeSequins sampBefore, sampAfter;
            
            std::map<ChrID, std::map<Locus, SampledInfo>> c2v;
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            Method meth = Method::Mean;
            
            // Defined only if meth==Prop
            Proportion p;
        };
        
        static Stats analyze(const FileName &gen,
                             const FileName &seqs,
                             const Options &o = Options());

        static void report(const std::vector<FileName> &,
                           const Options &o = Options());
    };
}

#endif