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
            Base trim = 2;
            
            // Defined only if meth==Prop
            Proportion p = NAN;
            
            // Defined only if meth==Reads
            Counts reads = NAN;

            // How calibration is calculated
            Method meth = Method::Mean;
        };
        
        enum class Status
        {
            ReverseReverse,
            ReverseNotMapped,
            ForwardForward,
            ForwardVarQuin,
            ForwardNotMapped,
            NotMappedNotMapped,
            ReverseHang,
            ForwardHang,
            ReverseLadQuin,
            LadQuinLadQuin,
            Original,
            Passed,
            Sampled
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
                counts[Status::ReverseReverse]  = 0;
                counts[Status::LadQuinLadQuin]  = 0;
                counts[Status::ReverseLadQuin]  = 0;
                counts[Status::ReverseHang] = 0;
                counts[Status::ForwardHang] = 0;
                counts[Status::ReverseNotMapped] = 0;
                counts[Status::ForwardForward] = 0;
                counts[Status::ForwardVarQuin] = 0;
                counts[Status::ForwardNotMapped] = 0;
                counts[Status::NotMappedNotMapped] = 0;
            }
            
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
            
            /*
             * Calibration statistics
             */
            
            struct Calibration
            {
                // Number of alignments for sequins (excluding LadQuin) before calibration
                Counts bSeqs = 0;
                
                // Number of alignments for sequins (excluding LadQuin) after calibration
                Counts aSeqs = 0;
                
                std::map<SequinID, Coverage> covs;
                std::map<SequinID, Proportion> norms;

                std::vector<double> allBeforeEndoC;
                std::vector<double> allBeforeSeqsC;
                
                // Required for summary statistics
                std::vector<double> allNorms;

                // Summary statistics for normalization factors
                inline double normMean() const
                {
                    return SS::mean(allNorms);
                }
                
                // Summary statistics for normalization factors
                inline double normSD() const
                {
                    return SS::getSD(allNorms);
                }

                // Average sequence coverage for endogenous before normalization
                inline double meanBEndo() const
                {
                    return SS::mean(allBeforeEndoC);
                }
                
                // Average sequence coverage for sequins before normalization
                inline double meanBSeqs() const
                {
                    return SS::mean(allBeforeSeqsC);
                }
            };
            
            Calibration cStats;
            
            /*
             * Mapping statistics
             */
            
            struct Mapping
            {
                // List of sequins
                std::set<SequinID> seqs;
                
                // Mapping between sequins to their endogenous regions
                std::map<SequinID, std::pair<ChrID, std::string>> s2e;
                
                // Mapping between sequins to their sequin regions
                std::map<SequinID, std::pair<ChrID, std::string>> s2s;
                
                // Alignment coverage for sequins (before and after) and endogenous regions
                ID2Intervals eInters, bInters, aInters;

                // Eg: GS_037 to chr2
                inline ChrID s2c(const SequinID &x) const
                {
                    return s2e.at(x).first;
                }
            };
            
            Mapping mStats;
            
            /*
             * General statistics
             */
            
            struct GeneralStats
            {
                // Number of reference regions
                Counts nRegs;
                
                // Total alignments before sampling for endogenous
                Counts bTEndo = 0;
                
                // Total alignments before sampling for sequins
                Counts bTSeqs = 0;
                
                // Total alignments after sampling for endogenous
                Counts aTEndo = 0;

                // Total alignments after sampling for sequins
                Counts aTSeqs = 0;

                // Regional alignments before sampling for endogenous
                Counts bREndo = 0;
                
                // Regional alignments after sampling for endogenous
                Counts aREndo = 0;
            };
            
            GeneralStats gStats;
            
            Coverage afterSeqs;
            
            // Alignment records for sequins (we'll need them for sampling)
            std::vector<ParserBAM::Data> s1, s2;

            // Sampling statistics
            std::map<SequinID, SampledInfo> c2v;
            
            std::map<VProcess::Status, Counts> counts;
        };

        static Stats analyze(const FileName &, const Options &);
        static void report(const FileName &, const Options &);
    };
}

#endif
