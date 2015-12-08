#ifndef COVERAGE_TOOL_HPP
#define COVERAGE_TOOL_HPP

#include "parsers/parser.hpp"
#include "stats/analyzer.hpp"
#include "data/intervals.hpp"
#include "data/alignment.hpp"

namespace Anaquin
{
    struct CoverageTool
    {
        struct Stats : public AlignmentStats, public SingleInputStats
        {
            Intervals<> inters;
        };
        
        struct Mapping
        {
            GenericID id;
            
            // Where the mapping occurs
            Locus l;
        };

        // Whether to proceed with the alignment
        typedef std::function<bool (const Alignment &, const ParserProgress &)> AlignFunctor;

        // Whether to proceed with the coverage
        typedef std::function<bool (const ChromoID &, Base, Base, Coverage)> CoverageFunctor;

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

            Interval::IntervalID id = Standard::chrT;
        };
        
        struct CoverageBedGraphOptions
        {
            // Where the bedgraph should be generated
            FileName file;
            
            // How the data should be written
            std::shared_ptr<Writer> writer;
        };
        
        // Analyze a BAM file sorted by position
        static Stats stats(const FileName &, AlignFunctor);

        // Analyze a list of mappings assume sorted
        static Stats stats(const std::vector<Mapping> &, const std::map<GenericID, Base> &);
        
        // Generate a bedgraph for the coverage
        static void bedGraph(const Stats &, const CoverageBedGraphOptions &, CoverageFunctor);

        // Generate summary statistics for the coverage
        static void summary(const Stats &, const CoverageReportOptions &, CoverageFunctor);
    };
}

#endif