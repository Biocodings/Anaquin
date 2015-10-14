#include "tools/coverage.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

CoverageTool::Stats CoverageTool::stats(const FileName &file, Functor f)
{
    CoverageTool::Stats stats;

    stats.src = file;
    
    /*
     * Reference: https://github.com/arq5x/bedtools2/blob/master/src/genomeCoverageBed/genomeCoverageBed.cpp
     */

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        auto addCoverage = [&](const ChromoID &id, Base start, Base end)
        {
            const auto size = stats.chroms.at(id).size;
            assert(size);
            
            // process the first line for this chromosome.
            // make sure the coordinates fit within the chrom
            if (start < size)
            {
                stats.chroms.at(id).covs[start].starts++;
            }
            if (end < size)
            {
                stats.chroms.at(id).covs[end].ends++;
            }
            else
            {
                stats.chroms.at(id).covs.back().ends++;
            }
        };
        
        if (!align.i)
        {
            if      (!align.mapped)                       { stats.unmapped++; }
            else if (align.id != Standard::instance().id) { stats.n_expT++;   }
            else                                          { stats.n_chrT++;   }
        }

        // Proceed with the alignment?
        if (f(align, info.p))
        {
            if (!stats.chroms.count(align.id))
            {
                stats.chroms[align.id].name = align.id;
                stats.chroms[align.id].size = info.size;
                stats.chroms[align.id].covs.resize(stats.chroms[align.id].size);
            }

            /*
             * This is like looping for blocks in AddBlockedCoverage()
             */

            addCoverage(align.id, align.l.start, align.l.end);
        }
    });

    // Find the depth for ith percentile
    auto percent = [&](Percentage p, const std::map<Base, Counts> &depths)
    {
        Percentage i = 0;
        
        for (const auto &depth : depths)
        {
            i += depth.second;
            
            // Have we reached our percentile?
            if (i >= p)
            {
                return depth.second;
            }
        }
        
        return depths.end()->second;
    };

    /*
     * Calculate descriptive statistics
     */
    
    for (auto &chrom : stats.chroms)
    {
        Counts n = 0;
        Coverage sums = 0;
        std::map<Base, Counts> depths;

        chrom.second.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            n++;
            sums += depth;
            depths[depth]++;
        });
        
        assert(n);
        
        stats.chroms[chrom.first].mean = sums / n;
        stats.chroms[chrom.first].min  = percent(0, depths);
        stats.chroms[chrom.first].max  = percent(n, depths);
        stats.chroms[chrom.first].p25  = percent(0.25 * n, depths);
        stats.chroms[chrom.first].p50  = percent(0.50 * n, depths);
        stats.chroms[chrom.first].p75  = percent(0.75 * n, depths);
    }
    
    return stats;
}

void CoverageTool::report(const CoverageTool::Stats &stats, const CoverageToolOptions &o)
{
    o.writer->open(o.bedGraph);

    // Coverage for chrT
    const auto &chrT = stats.chroms.at(Standard::instance().id);

    chrT.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
    {
        o.writer->write((boost::format("%1%\t%2%\t%3%\t%4%") % id
                                                             % i
                                                             % j
                                                             % depth).str());
    });
    
    o.writer->close();

    /*
     * Generating summary statistics
     */

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Experiment: %2%\n"
                         "   Synthetic: %3%\n\n"
                         "   Reference: %4%\n"
                         "   Reference Bases: %5%\n\n"
                         "   Minimum: %6%\n"
                         "   Maximum: %7%\n"
                         "   Mean:    %8%\n"
                         "   25th: %9%\n"
                         "   50th: %10%\n"
                         "   75th: %11%\n"
    ;

    o.writer->open(o.summary);
    o.writer->write((boost::format(summary) % stats.src
                                            % stats.n_expT
                                            % stats.n_chrT
                                            % o.refs
                                            % o.size
                                            % chrT.min
                                            % chrT.max
                                            % chrT.mean
                                            % chrT.p25
                                            % chrT.p50
                                            % chrT.p75).str());
    o.writer->close();
}