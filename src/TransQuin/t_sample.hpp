#ifndef T_SAMPLE_HPP
#define T_SAMPLE_HPP

#include "tools/sample.hpp"
#include "tools/coverage.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TSample
    {
        struct SampleImpl : public Subsampler::SampleImpl
        {
            #define V_REF Standard::instance().r_var
            
            inline ChrID genoID() const
            {
                return V_REF.genoID();
            }
            
            inline bool shouldGenomic(const ChrID &id, const Locus &l) const
            {
                return static_cast<bool>(V_REF.findGeno(id, l));
            }
            
            inline bool shouldChrT(const ChrID &id, const Locus &l) const
            {
                return static_cast<bool>(V_REF.match(l, MatchRule::Contains));
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
            const auto stats = VSample::stats(file, o);
            //Subsampler::report(file, o);
        }
    };
}

#endif