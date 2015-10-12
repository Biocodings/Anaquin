#include "analyzers/coverage.hpp"

using namespace Anaquin;

CoverageTool::Stats CoverageAnalyzer::stats(const FileName &file, const Options &o)
{
    const auto stats = CoverageTool::analyze(file, [&](const Alignment &align, const ParserProgress &p)
    {
        return true;
    });
    
    return stats;
}

void CoverageAnalyzer::report(const FileName &file, const FileName &output, const Options &o)
{
    const auto stats = CoverageAnalyzer::stats(file, o);

    Base depth = 0;
    long lastStart = -1;
    long lastDepth = -1;

    o.writer->open(output);
    
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