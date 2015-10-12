#ifndef COVERAGE_HPP
#define COVERAGE_HPP

#include <vector>
#include "parsers/parser.hpp"
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    struct CoverageTool
    {
        struct Depth
        {
            Base starts;
            Base ends;
        };

        struct Chromosome
        {
            ChromoID name;

            // Size of the chromosome
            Base size;

            // Base coverage for the chromosome
            std::vector<Depth> covs;
        };
        
        struct Stats
        {
            std::map<ChromoID, Chromosome> chroms;
        };

        // Whether to proceed with the alignment
        typedef std::function<bool (const Alignment &, const ParserProgress &)> Functor;

        // Analyze a BAM file sorted by position
        static Stats analyze(const FileName &, Functor);

        // Report a BAM file sorted by position
        static void report(const FileName &, const FileName &, Functor);
    };
}

#endif