#ifndef GI_M_ABUNDANCE_HPP
#define GI_M_ABUNDANCE_HPP

namespace Spike
{
    struct MAssemblyStats
    {
        // Percentage of the genome that is contained in the assembly
        Percentage cov_genome;

        // Percentage of the genes in the genome that are contained in the assembly
        std::map<SequinID, Percentage> cov_gene;

        // Empirical distribution for contigs  
        Data contigs;
    };

    struct MAssembly
    {
        struct Options : public AnalyzerOptions
        {
            // Empty Implementation
        };

        static MAssemblyStats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif