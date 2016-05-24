#ifndef M_SAMPLE_HPP
#define M_SAMPLE_HPP

#include "tools/sample.hpp"
#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MSample
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

            inline bool shouldSynthetic(const ChrID &id, const Locus &l) const
            {
                return static_cast<bool>(REF.match(l, MatchRule::Contains));
            }
        };
        
        struct ReportImpl : public Subsampler::ReportImpl
        {
            inline FileName summary() const
            {
                return "MetaSubsample_summary.stats";
            }
            
            virtual FileName beforeBG() const
            {
                return "MetaSubsample_before.bedgraph";
            }
            
            virtual FileName afterBG() const
            {
                return "MetaSubsample_after.bedgraph";
            }
            
            virtual FileName sampled() const
            {
                return "MetaSubsample_sampled.sam";
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