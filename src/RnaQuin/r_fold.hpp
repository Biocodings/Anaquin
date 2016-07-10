/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_FOLD_HPP
#define R_FOLD_HPP

#include "data/dtest.hpp"
#include <boost/format.hpp>
#include "stats/analyzer.hpp"

// Defined in resources.cpp
extern Anaquin::Scripts PlotTLODR();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTROC();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTMA();

// Defined in resources.cpp
extern Anaquin::FileName MixRef();

// Defined in resources.cpp
extern Anaquin::FileName GTFRef();

namespace Anaquin
{
    struct RFold : public Analyzer
    {
        template <typename Options> static void generateLODR(const FileName &file,
                                                             const FileName &csv,
                                                             const Options &o)
        {
            o.generate(file);
            o.writer->open(file);
            o.writer->write(RWriter::createScript(csv, PlotTLODR()));
            o.writer->close();
        }

        enum class Metrics
        {
            Gene,
            Isoform
        };

        struct Options : public DoubleMixtureOptions
        {
            Options() {}
            Metrics metrs;
        };

        struct Stats : public LinearStats, public MappingStats, public AnalyzerStats
        {
            struct Data
            {
                // Expcted log-fold ratio
                Concent exp;
                
                // Measured log-fold ratio
                Concent obs;
                
                // Standard deviation
                double se;
                
                // Base mean
                double mean;
                
                Probability p, q;
            };
            
            std::map<SequinID, Data> data;
            
            /*
             * Optional inputs
             */
            
            // Normalized average counts for the replicates
            std::vector<double> means;
            
            // Log-fold ratios standard deviation
            std::vector<double> ses;
            
            // Average counts for each condition if provided
            //std::vector<std::map<std::string, Counts>> avgs;
        };

        /*
         * Classify and construct a vector of TP/FP/TN/FN, given the q-values and expected
         * fold-changes. The threshold for the TP classification is also required.
         */
        
        static std::vector<std::string> classify(const std::vector<double> &,
                                                 const std::vector<double> &,
                                                 double qCut,
                                                 double foldCut);

        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif