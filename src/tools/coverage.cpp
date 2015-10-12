#include "tools/coverage.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

void CoverageTool::report(const FileName &file, const FileName &bg, Functor f)
{
    const auto stats = CoverageTool::analyze(file, f);

    Base depth = 0;
    long lastStart = -1;
    long lastDepth = -1;

    FileWriter writer("");
    writer.open(file);

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
                    writer.write((boost::format("%1%\t%2%\t%3%\t%4%") % chr.name
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
             * changes. Update
             */
            
            // Default: the depth has not changed, so we will not print anything.
            // Proceed until the depth changes.
            // Update depth
            depth = depth - chr.covs[j].ends;
        }
        
        // Print information about the last position
        if ((lastDepth != -1) && (lastDepth > 0))
        {
            writer.write((boost::format("%1%\t%2%\t%3%\t%4%") % chr.name
                                                              % lastStart
                                                              % lastDepth).str());
            //std::cout << chr.name << "\t" << lastStart << "\t" << chr.size << "\t" << lastDepth << std::endl;
        }
    }
    
    writer.close();
}

CoverageTool::Stats CoverageTool::analyze(const FileName &file, Functor f)
{
    CoverageTool::Stats stats;

    /*
     * Reference: https://github.com/arq5x/bedtools2/blob/master/src/genomeCoverageBed/genomeCoverageBed.cpp
     */
    
    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &p)
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
                assert(false);
                //_currChromCoverage[_currChromSize-1].ends++;
            }
        };
        
        // Proceed with the alignment?
        if (f(align, p))
        {
            if (!stats.chroms.count(align.id))
            {
                stats.chroms[align.id].name = align.id;
                stats.chroms[align.id].size = 8457082;
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