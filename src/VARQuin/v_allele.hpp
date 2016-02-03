#ifndef V_ALLELE_HPP
#define V_ALLELE_HPP

#include "stats/analyzer.hpp"
#include "VARQuin/VARQuin.hpp"

namespace Anaquin
{
    struct VAllele
    {
        struct Options : public AnalyzerOptions
        {
            Caller caller;
        };
        
        struct Stats
        {
            typedef LinearStats Data;

            std::map<ChromoID, Data> data;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif