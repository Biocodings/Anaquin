#ifndef V_PROCESS_HPP
#define V_PROCESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VProcess
    {
        enum class Method
        {
            Mean,
            Median,
            Reads,
            Prop,
        };

        struct Options : AnalyzerOptions
        {
            // How much edge effects?
            Base trim = 1;
            
            // Defined only if meth==Prop
            Proportion p = NAN;
            
            // Defined only if meth==Reads
            Counts reads = NAN;

            Method meth = Method::Mean;
        };
        
        enum class Status
        {
            ReverseReverse,
            ReverseNotMapped,
            ForwardForward,
            ForwardReverse,
            ForwardNotMapped,
            NotMappedNotMapped,
            RevHang,
            ForHang
        };

        struct SampledInfo
        {
            SequinID rID;
            
            // Alignment coverage for the endogenous sample
            Coverage endo;
            
            // Alignment coverage before subsampling
            Coverage before;
            
            // Alignment coverage after subsampling
            Coverage after;
            
            // Number of alignments before and after
            Counts nEndo, nBefore, nAfter;
            
            // Normalization factor
            Proportion norm;
        };

        struct Stats : public MappingStats
        {
            std::vector<double> allBeforeEndoC;
            std::vector<double> allBeforeSeqsC;

            // Required for summary statistics
            std::vector<double> allNorms;

            // Sampling statistics
            std::map<ChrID, std::map<Locus, SampledInfo>> c2v;
            
            // Alignment coverage for sequins and genome
            ID2Intervals sInters, gInters;
            
            std::map<VProcess::Status, Counts> counts;
        };

        static Stats analyze(const FileName &, const Options &);
        static void report(const FileName &, const Options &);
    };
}

#endif
