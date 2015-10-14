#ifndef COVERAGE_TOOL_HPP
#define COVERAGE_TOOL_HPP

#include <vector>
#include "parsers/parser.hpp"
#include "stats/analyzer.hpp"
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

        struct ChromoCoverage
        {
            ChromoID name;

            // Size of the chromosome
            Base size;

            // Base coverage for the chromosome
            std::vector<Depth> covs;
            
            // Minimum coverage
            Coverage min;
            
            // Maximum coverage
            Coverage max;
            
            // Average coverage
            Coverage mean;
            
            // 25th percentile coverage
            Coverage p25;
            
            // 50th percentile coverage
            Coverage p50;
            
            // 75th percentile coverage
            Coverage p75;
            
            /*
             * Define a template function for looping over the coverage for each chromosome
             * in BedGraph format.
             */

            template <typename T> void bedGraph(T t) const
            {
                Base depth = 0;
                long lastStart = -1;
                long lastDepth = -1;

                for (auto j = 0; j < size; j++)
                {
                    depth += covs[j].starts;

                    if (depth != lastDepth)
                    {
                        /*
                         * Coverage depth has changed, print the last interval coverage (if any)
                         *
                         * Print if:
                         *
                         *   (1) depth>0  (the default running mode),
                         *   (2) depth==0 and the user requested to print zero covered regions
                         */
                        
                        if ((lastDepth != -1) && (lastDepth > 0))
                        {
                            t(name, lastStart, j, lastDepth);
                        }
                        
                        // Set current position as the new interval start + depth
                        lastStart = j;
                        lastDepth = depth;
                    }
                    
                    depth = depth - covs[j].ends;
                }

                // Print information about the last position
                if ((lastDepth != -1) && (lastDepth > 0))
                {
                    t(name, lastStart, size, lastDepth);
                }
            }
        };

        struct Stats : public AlignmentStats
        {
            std::map<ChromoID, ChromoCoverage> chroms;
        };

        // Whether to proceed with the alignment
        typedef std::function<bool (const Alignment &, const ParserProgress &)> Functor;

        struct CoverageToolOptions
        {
            // Filename for the summary statistics
            FileName summary;
            
            // Filename for the generated bedgraph
            FileName bedGraph;
            
            // Number of sequins
            Counts refs;
            
            // Size of all sequins
            Base size;
            
            // Where the data should be written
            std::shared_ptr<Writer> writer;
        };
        
        // Analyze a BAM file sorted by position
        static Stats stats(const FileName &, Functor);

        // Report a BAM file sorted by position
        static void report(const Stats &, const CoverageToolOptions &);
    };
}

#endif