#ifndef V_SAMPLE_HPP
#define V_SAMPLE_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_bambed.hpp"

namespace Anaquin
{
    typedef std::map<ChrID, std::map<Locus, Proportion>> NormFactors;
    
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
        
        struct CalibrateStats
        {
            Counts nEndo = 0;
            Counts nSeqs = 0;
         
            ParserBAMBED::Stats es;
            ParserBAMBED::Stats ss;
            
            std::vector<double> allBeforeEndoC;
            std::vector<double> allBeforeSeqsC;
            
            // Required for summary statistics
            std::vector<double> allNorms;

            // Normalization for each region
            NormFactors norms;

            std::map<ChrID, std::map<Locus, SampledInfo>> c2v;
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

            // Defined only if meth==Prop
            Proportion p = NAN;
            
            // Defined only if meth==Reads
            Counts reads = NAN;
        };
        
        static ParserBAMBED::Stats sample(const FileName    &,
                                          const NormFactors &,
                                          VSample::Stats    &,
                                          const C2Intervals &,
                                          const C2Intervals &,
                                          const VSample::Options &);
        
        static CalibrateStats check(const FileName &,
                                    const FileName &,
                                    const C2Intervals &,
                                    const C2Intervals &,
                                    const VSample::Options &);

        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
