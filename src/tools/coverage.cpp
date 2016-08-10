#include "tools/coverage.hpp"
#include "writers/file_writer.hpp"

// Defined in main.cpp
extern bool __showInfo__;

using namespace Anaquin;

CoverageTool::Stats CoverageTool::stats(const FileName &file, std::map<ChrID, Intervals<>> &inters)
{
    CoverageTool::Stats stats;
 
    stats.src = file;
    
    ParserSAM::parse(file, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (__showInfo__ && info.p.i && !(info.p.i % 1000000))
        {
            std::cout << std::to_string(info.p.i) << std::endl;
        }
        
        if (Standard::isSynthetic(x.cID))
        {
            stats.countSyn++;
        }
        else if (x.cID != "*")
        {
            stats.countGen++;
        }
        else
        {
            stats.countNA++;
        }
        
        if (x.mapped && inters.count(x.cID))
        {
            if (info.skip)
            {
                return;
            }
            
            Locus l;
            bool spliced;
            bool added = false;
            
            while (x.nextCigar(l, spliced))
            {
                std::vector<Interval *> mats;
                if (inters[x.cID].overlap(l, &mats))
                {
                    for (const auto &i : mats)
                    {
                        i->map(l);
                        
                        if (!added)
                        {
                            added = true;
                            stats.hist[x.cID]++;
                        }
                    }
                }
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
            o.writer->write((boost::format("%1%\t%2%\t%3%\t%4%") % id % i % j % depth).str());
        }
    });

    o.writer->close();
}