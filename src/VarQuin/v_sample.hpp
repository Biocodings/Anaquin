#ifndef V_SAMPLE_HPP
#define V_SAMPLE_HPP

#include "tools/sample.hpp"
#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VSample
    {
        struct SampleImpl : public Subsampler::StatsImpl
        {
            #define REF Standard::instance().r_var            
        };
        
        struct ReportImpl : public Subsampler::ReportImpl
        {
            inline FileName summary() const
            {
                return "VarSubsample_summary.stats";
            }
            
            virtual FileName beforeBG() const
            {
                return "VarSubsample_before.bedgraph";
            }
            
            virtual FileName afterBG() const
            {
                return "VarSubsample_after.bedgraph";
            }
            
            virtual FileName sampled() const
            {
                return "VarSubsample_sampled.sam";
            }
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}

            Subsampler::CoverageMethod method = Subsampler::CoverageMethod::ArithAverage;
        };

        typedef Subsampler::Stats Stats;

        static Stats stats(const FileName &file, const Options &o = Options())
        {
            return Subsampler::stats(file, o, SampleImpl());
        }

        static void report(const FileName &file, const Options &o = Options())
        {
            Subsampler::report(file, o, SampleImpl(), ReportImpl());
        }
    };
}

#endif