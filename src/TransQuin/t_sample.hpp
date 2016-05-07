#ifndef T_SAMPLE_HPP
#define T_SAMPLE_HPP

#include "tools/sample.hpp"
#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TSample
    {
        struct SampleImpl : public Subsampler::StatsImpl
        {
            #define REF Standard::instance().r_var
            
            inline ChrID genoID() const
            {
                return REF.genoID();
            }
            
            inline bool shouldGenomic(const ChrID &id, const Locus &l) const
            {
                return static_cast<bool>(REF.findGeno(id, l));
            }
            
            inline bool shouldChrT(const ChrID &id, const Locus &l) const
            {
                return static_cast<bool>(REF.match(l, MatchRule::Contains));
            }
        };
        
        struct ReportImpl : public Subsampler::ReportImpl
        {
            inline FileName summary() const
            {
                return "TransSubsample_summary.stats";
            }
            
            virtual FileName beforeBG() const
            {
                return "TransSubsample_before.bedgraph";
            }
            
            virtual FileName afterBG() const
            {
                return "TransSubsample_after.bedgraph";
            }
            
            virtual FileName sampled() const
            {
                return "TransSubsample_sampled.sam";
            }
            
            virtual Counts countSeqs() const
            {
                return REF.countSeqs();
            }
            
            virtual Counts countInters() const
            {
                return REF.countInters();
            }
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            // How coverage is calculated
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