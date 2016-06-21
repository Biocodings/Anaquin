#include "tools/coverage.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

CoverageTool::Stats CoverageTool::stats__(const FileName &file, std::map<ChrID, Intervals<>> &inters)
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

CoverageTool::Stats CoverageTool::stats(const FileName &file, AlignFunctor f)
{
    CoverageTool::Stats stats;

    stats.src  = file;
    
    ParserSAM::parse(file, [&](const Alignment &x, const ParserSAM::Info &info)
    {
        stats.update(x);

        if (x.mapped && f(x, info.p))
        {
            if (!stats.inters.find(x.cID))
            {
                // Add a new interval for the chromosome (giving the length of the chromosome)
                stats.inters.add(Interval(x.cID, Locus(0, info.length-1)));
            }

            stats.hist[x.cID]++;
            stats.inters.find(x.cID)->add(x.l);
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

//void CoverageTool::bedGraph(const Stats &stats, const CoverageBedGraphOptions &o, CoverageFunctor f)
//{
//    o.writer->open(o.file);
//
//    for (const auto &i : stats.inters.data())
//    {
//        const auto chr = i.second;
//
//        chr.bedGraph([&](const ChrID &id, Base i, Base j, Base depth)
//        {
//            if (depth && f(id, i, j, depth))
//            {
//                o.writer->write((boost::format("%1%\t%2%\t%3%\t%4%") % id
//                                                                     % i
//                                                                     % j
//                                                                     % depth).str());
//            }
//        });
//    }
//
//    o.writer->close();
//}