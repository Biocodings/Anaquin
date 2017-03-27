#ifndef V_VREPORT_HPP
#define V_VREPORT_HPP

#include "tools/system.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VVReport : public Analyzer
    {
        struct Options : public AnalyzerOptions
        {
            FileName rMix, rVCF, rBED;
        };

        static void report(const FileName &hBAM, const FileName &sBAM, const FileName &qVCF, const Options &o)
        {
            const auto code = ((boost::format("library(Anaquin); createVReport('%1%', '%2%', '%3%', '%4%', '%5%', '%6%')"))
                                      % o.rMix
                                      % o.rVCF
                                      % o.rBED
                                      % hBAM
                                      % sBAM
                                      % qVCF).str();
            System::runRScript(code);
        }
    };
}

#endif
