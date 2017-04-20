#ifndef V_SAMPLE_HPP
#define V_SAMPLE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VSample
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
            
            // Alignment coverage for the endogenous sample
            Coverage endo;
            
            // Alignment coverage before subsampling
            Coverage before;
            
            // Alignment coverage after subsampling
            Coverage after;
            
            // Normalization factor
            Proportion norm;
        };
        
        struct GenomeSequins
        {
            Counts nEndo = 0;
            Counts nSeqs = 0;
        };
        
        struct Stats
        {
            // Total number of subsampling regions
            Counts count = 0;
            
            // Number of regions without alignment (genomic)
            Counts noGAlign = 0;
            
            // Number of regions without alignment (synthetic)
            Counts noSAlign = 0;
            
            Coverage afterEndo,  afterSeqs;
            Coverage beforeEndo, beforeSeqs;
            
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
            
            Base edge = 0;

            // Trim alignment edge effects?
            Base trim = 0;
            
            // Defined only if meth==Prop
            Proportion p = NAN;
            
            // Defined only if meth==Reads
            Counts reads = NAN;
        };
        
        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
