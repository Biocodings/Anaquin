#include "tools/coverage.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

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