#ifndef COVERAGE_TOOL_HPP
#define COVERAGE_TOOL_HPP

#include "data/intervals.hpp"
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    struct CoverageTool
    {
        struct Stats : public AlignmentStats, public SequinStats
        {
            FileName src;
        };

        struct CoverageReportOptions
        {
            // Filename for the summary statistics
            FileName summary;
            
            // Number of sequins
            Counts refs;

            // Size of all sequins
            Base length;
            
            // Where the data should be written
            std::shared_ptr<Writer> writer;

            Interval::IntervalID id = ChrIS;
        };
        
        struct CoverageBedGraphOptions
        {
            // Where the bedgraph should be generated
            FileName file;
            
            // How the data should be written
            std::shared_ptr<Writer> writer;
        };
        
        static Stats stats(const FileName &, std::map<ChrID, Intervals<>> &);

        static void bedGraph(const ID2Intervals &, const CoverageBedGraphOptions &);
    };
}

#endif