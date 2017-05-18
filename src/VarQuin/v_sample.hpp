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
            SequinID rID;
            
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

        struct GenomeSequins
        {
            Counts nEndo = 0;
            Counts nSeqs = 0;
        };
        
        struct Stats
        {
            CalibrateStats cStats;
            
            Coverage afterSeqs;
            
            GenomeSequins tBefore, tAfter;
            GenomeSequins sBefore, sAfter;
            
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

        static GenomeSequins tBefore(const CalibrateStats &, const ParserBAMBED::Stats &);
        static GenomeSequins tAfter (const CalibrateStats &, const ParserBAMBED::Stats &);
        static GenomeSequins sBefore(const CalibrateStats &, const ParserBAMBED::Stats &);
        static GenomeSequins sAfter (const CalibrateStats &, const ParserBAMBED::Stats &);
        
        static double afterSeqsC(const C2Intervals &tRegs, std::map<ChrID, std::map<Locus, SampledInfo>> &c2v, VSample::Options o);
        
        static ParserBAMBED::Stats sample(const FileName    &,
                                          const NormFactors &,
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
    
    inline std::string meth2Str(VSample::Method meth)
    {
        switch (meth)
        {
            case VSample::Method::Mean:   { return "Mean";       }
            case VSample::Method::Median: { return "Median";     }
            case VSample::Method::Reads:  { return "Reads";      }
            case VSample::Method::Prop:   { return "Proportion"; }
        }
    };
}

#endif
