#ifndef M_ABUNDANCE_HPP
#define M_ABUNDANCE_HPP

#include "meta/m_assembly.hpp"

namespace Anaquin
{
    struct MAbundance
    {
        template <typename Options> static void logOptions(const Options &o)
        {
            switch (o.coverage)
            {
                case WendyAlgorithm: { o.info("Wendy Smoothing");           break; }
                case KMerCov_Contig: { o.info("K-mer coverage per contig"); break; }
                case KMerCov_Sequin: { o.info("K-mer coverage per sequin"); break; }
            }
        };

        struct Stats : public LinearStats, public MappingStats
        {
            Sensitivity ss;

            // Distribution of the sequins
            SequinHist h = Standard::instance().r_meta.hist();
        };

        enum CoverageMethod
        {
            WendyAlgorithm,            
            KMerCov_Contig, // K-mer coverage relative to the size of the contig
            KMerCov_Sequin, // K-mer coverage relative to the size of the sequin
        };

        struct Options : public MAssembly::Options
        {
            // Required by the GCC compiler ...
            Options() {}

            CoverageMethod coverage = WendyAlgorithm;
        };
        
        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif