#ifndef GI_PARSER_HPP
#define GI_PARSER_HPP

#include <memory>
#include "writers/writer.hpp"
#include "exceptions/invalid_extension.hpp"

namespace Spike
{
    struct ParserStats
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
    
    struct ParserOptions
    {
        std::shared_ptr<Writer> writer;
    };
}

#endif