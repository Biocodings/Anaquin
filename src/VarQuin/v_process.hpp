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
            Stats()
            {
                counts[Status::RevHang] = 0;
                counts[Status::ForHang] = 0;
                counts[Status::LadQuin] = 0;
                counts[Status::ReverseReverse] = 0;
                counts[Status::ForwardForward] = 0;
                counts[Status::ForwardReverse] = 0;
                counts[Status::ReverseLadQuin] = 0;
                counts[Status::ReverseNotMapped] = 0;
                counts[Status::ForwardNotMapped] = 0;
                counts[Status::NotMappedNotMapped] = 0;
            }
            
            /*
             * General statistics
             */
            
            // Number of reference regions
            Counts nRegs;

            struct Trim
            {
                // Number of reads trimmed on left
                Counts left = 0;
                
                // Number of reads trimmed on right
                Counts right = 0;
                
                // Number of alignments before trimming
                Counts before = 0;
                
                // Number of alignments after trimming
                Counts after  = 0;
            };
            
            struct Ladder
            {
                // Number of alignments for ladders
                Counts nLad = 0;
            };
            
            // Trimming statistics
            Trim trim;
            
            // Ladder statistics
            Ladder lad;
            
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
