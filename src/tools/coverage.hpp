#ifndef COVERAGE_TOOL_HPP
#define COVERAGE_TOOL_HPP

#include "parsers/parser.hpp"
#include "stats/analyzer.hpp"
#include "data/intervals.hpp"
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    struct CoverageTool
    {
        struct Stats : public AlignmentStats, public SequinStats
        {
            FileName src;

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
        typedef std::function<bool (const ChrID &, Base, Base, Coverage)> CoverageFunctor;

        struct CoverageReportOptions
        {
            // Filename for the summary statistics
            FileName summary;
            
            FileName rAnnot, rGeno;
            
            // Number of sequins
            Counts refs;

            // Size of all sequins
            Base length;
            
            // Where the data should be written
            std::shared_ptr<Writer> writer;

            Interval::IntervalID id = ChrT;
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

        static Stats stats__(const FileName &, std::map<ChrID, Intervals<>> &);

        // Analyze a BAM file sorted by position
        template <typename F> static Stats stats_(const FileName &file, const Hist &hist, F f)
        {
            CoverageTool::Stats stats;
            
            stats.src  = file;
            stats.hist = hist;
            
            ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::Info &info)
            {
                stats.update(align);
                
                if (align.mapped)
                {
                    const auto match = f(align, info.p);
                    
                    if (match)
                    {
                        if (!stats.inters.find(align.cID))
                        {
                            // Add a new interval for the chromosome
                            stats.inters.add(Interval(align.cID, Locus(0, info.length-1)));
                        }
                        
                        //stats.hist.at(match->id)++;
                        stats.inters.find(align.cID)->add(align.l);
                    }
                }
            });
            
            return stats;
        }

        // Analyze a list of mappings assume sorted
        static Stats stats(const std::vector<Mapping> &, const std::map<GenericID, Base> &);
        
        // Generate a bedgraph for the coverage
        static void bedGraph(const Stats &, const CoverageBedGraphOptions &, CoverageFunctor);

        // Generate summary statistics for the coverage
        static void summary(const Stats &, const CoverageReportOptions &, CoverageFunctor);

        static Scripts writeCSV(const Stats &, const CoverageReportOptions &);
    };
}

#endif