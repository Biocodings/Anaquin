#ifndef GI_F_EXPRESS_HPP
#define GI_F_EXPRESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    /*
     * This analyzer is designed to analyse the relative abundance between fusion genes
     * and normal genes. For example, if we have a fusion gene with known abundance of
     * 100 and a normal gene with known abundance of 10. The expected relative ratio would
     * be 100/10 = 10. How does it compare with the measured relative ratio?
     *
     * The inputs will be a genes tracking file from cufflinks. The file will give the measured
     * RPKM for each gene.
     */

    struct FExpress
    {
        struct Stats : public ModelStats
        {
            // Empty Implementation            
        };
        
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif