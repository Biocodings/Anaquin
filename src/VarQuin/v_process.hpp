#ifndef V_PROCESS_HPP
#define V_PROCESS_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_bam.hpp"

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
            Base edge = 0;
            
            // Should we trim sequin reads?
            bool shouldTrim = true;
            
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
            ForHang,
            ReverseLadQuin,
            LadQuin
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
            // Number of alignments for ladders
            Counts nLad = 0;
            
            // Alignment records for sequins (we'll need them for sampling)
            std::vector<ParserBAM::Data> s1, s2;

            std::vector<double> allBeforeEndoC;
            std::vector<double> allBeforeSeqsC;

            // Required for summary statistics
            std::vector<double> allNorms;

            typedef std::map<ChrID, std::map<Locus, Proportion>> NormFactors;
            
            // Normalization for each region
            NormFactors norms;

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
