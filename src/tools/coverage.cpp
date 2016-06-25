#include "tools/coverage.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

CoverageTool::Stats CoverageTool::stats(const FileName &file, std::map<ChrID, Intervals<>> &inters)
{
    CoverageTool::Stats stats;
 
    stats.src = file;
    
    ParserSAM::parse(file, [&](const Alignment &x, const ParserSAM::Info &info)
    {
        stats.update(x);

        if (x.mapped && inters.count(x.cID))
        {
            const auto m = inters[x.cID].contains(x.l);
            
            if (m)
            {
                m->map(x.l);
                stats.hist[x.cID]++;
            }
        }
    });

    return stats;
}

void CoverageTool::bedGraph(const ID2Intervals &inters, const CoverageBedGraphOptions &o)
{
    o.writer->open(o.file);
    
    inters.bedGraph([&](const ChrID &id, Base i, Base j, Base depth)
    {
        if (depth)
        {
            o.writer->write((boost::format("%1%\t%2%\t%3%\t%4%") % id
                                                                 % i
                                                                 % j
                                                                 % depth).str());
        }
    });

    o.writer->close();
}