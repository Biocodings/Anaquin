#ifndef GI_ANALYZER_HPP
#define GI_ANALYZER_HPP

#include <memory>
#include "types.hpp"
#include "writers/writer.hpp"

namespace Spike
{
    struct AnalyzerStats
    {
        // Percentage of reads aligned with the reference chromosome
        Percentage pr;
        
        // Percentage of reads aligned with the query samples
        Percentage pq;
        
        // Total number of reads aligned
        Reads n = 0;
        
        // Number of reads aligned to the chromosome
        Reads nr = 0;
        
        // Number of reads aligned to the real sample
        Reads nq = 0;
        
        inline Percentage dilution() const
        {
            return nq ? static_cast<Percentage>(nr / nq) : 1;
        }
    };
    
    struct AnalyzerOptions
    {
        std::shared_ptr<Writer> writer;
    };
}

#endif