#ifndef V_DISCOVER_HPP
#define V_DISCOVER_HPP

#include "stats/analyzer.hpp"
#include "VARQuin/VARQuin.hpp"

namespace Anaquin
{
    struct VDiscover
    {
        struct Options : public AnalyzerOptions
        {
            Caller caller;
        };
        
        struct Stats
        {
            struct Classifed
            {
                // Expected allele frequency
                double eAllFreq;
                
                // Expected fold-change
                double eFold;
                
                CalledVariant query;

                // The sequin where it occurs
                const Variant *seq;
            };

            struct Data : public MappingStats
            {
                std::vector<Classifed> fps;
                std::vector<Classifed> tps;
                
                //long detected = 0;
                
                // Overall performance
                //Confusion m;
                
                //SequinHist h = Standard::instance().r_var.hist();
            };
            
            std::map<ChromoID, Data> data;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif