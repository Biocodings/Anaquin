#ifndef M_ABUNDANCE_HPP
#define M_ABUNDANCE_HPP

#include "meta/m_assembly.hpp"

namespace Anaquin
{
    struct MAbundance
    {
        struct Stats : public LinearStats, public MappingStats
        {
            Sensitivity ss;

            // Distribution of the sequins
            SequinHist h = Standard::instance().r_meta.hist();
        };

        enum CoverageMethod
        {
            KMerCov_Contig, // K-mer coverage relative to the size of the contig
            KMerCov_Sequin, // K-mer coverage relative to the size of the sequin
        };

        struct Options : public MAssembly::Options
        {
            // Required by the GCC compiler ...
            Options() {}

            CoverageMethod coverage = KMerCov_Contig;
        };
        
        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif