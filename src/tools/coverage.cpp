#include "tools/coverage.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

CoverageTool::Stats CoverageTool::stats(const FileName &file, Functor f)
{
    CoverageTool::Stats stats;

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
            else if (align.id != Standard::instance().id) { stats.n_hg38++;   }
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

    return stats;
}

void CoverageTool::report(const CoverageTool::Stats &stats, const CoverageToolOptions &o)
{
    Base depth = 0;
    long lastStart = -1;
    long lastDepth = -1;
    
    /*
     * Generating bedgraph for the standards
     */
    
    o.writer->open(o.bedGraph);
    
    for (auto i : stats.chroms)
    {
        const auto &chr = i.second;
        
        for (auto j = 0; j < chr.size; j++)
        {
            depth += chr.covs[j].starts;
            
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
                    o.writer->write((boost::format("%1%\t%2%\t%3%\t%4%") % chr.name
                                     % lastStart
                                     % j
                                     % lastDepth).str());
                    //std::cout << chr.name << "\t" << lastStart << "\t" << j << "\t" << lastDepth << std::endl;
                }
                
                // Set current position as the new interval start + depth
                lastDepth = depth;
                lastStart = j;
            }
            
            /*
             * Default: the depth has not changed, so we will not print anything. Proceed until the depth
             * changes.
             */
            
            depth = depth - chr.covs[j].ends;
        }
        
        // Print information about the last position
        if ((lastDepth != -1) && (lastDepth > 0))
        {
            o.writer->write((boost::format("%1%\t%2%\t%3%\t%4%") % chr.name
                             % lastStart
                             % lastDepth).str());
            //std::cout << chr.name << "\t" << lastStart << "\t" << chr.size << "\t" << lastDepth << std::endl;
        }
    }
    
    o.writer->close();
}