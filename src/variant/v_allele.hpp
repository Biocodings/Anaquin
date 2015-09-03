#ifndef GI_V_ALLELE_HPP
#define GI_V_ALLELE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAllele
    {
        typedef SingleMixtureOptions Options;
        typedef std::map<Variation, Counts> VarHist;

        struct Stats : public LinearStats
        {
            // Overall performance
            Confusion m;

            // The proportion of variations with alignment coverage
            double covered;

            // Measure of variant detection independent to sequencing depth or coverage
            double efficiency;

            // Distribution for the variants
            VarHist h = Analyzer::hist(Standard::instance().__v_vars__);
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif