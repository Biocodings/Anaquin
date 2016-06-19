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

            inline bool shouldGenomic(const ChrID &id, const Locus &l) const
            {
                return static_cast<bool>(REF.findGeno(id, l));
            }

            inline bool shouldSynthetic(const ChrID &id, const Locus &l) const
            {
                return static_cast<bool>(REF.match(l, MatchRule::Contains));
            }
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

            // How coverage is calculated. The default method follows Ira Deveson's 2016 paper.
            Subsampler::CoverageMethod method = Subsampler::CoverageMethod::Median;
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