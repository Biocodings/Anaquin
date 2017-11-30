#include <set>
#include <fstream>
#include "data/locus.hpp"
#include "tools/system.hpp"
#include "tools/bedtools.hpp"
#include "parsers/parser_bed.hpp"

using namespace Anaquin;

FileName BedTools::intersect(const FileName &x, const FileName &y, Base edge)
{
    std::set<ChrID> cs;
    std::map<ChrID, std::vector<ParserBed::Data>> m1;
    std::map<ChrID, std::vector<ParserBed::Data>> m2;
    
    ParserBed::parse(x, [&](ParserBed::Data &i, const ParserProgress &)
    {
        i.l.start += edge;
        i.l.end   -= edge;
        
        if (i.l.start > i.l.end)
        {
            throw std::runtime_error(i.name + " has length " + std::to_string(i.l.length()) + " , but the edge paramter is " + std::to_string(edge));
        }
        
        cs.insert(i.cID);
        m1[i.cID].push_back(i);
    });

    ParserBed::parse(y, [&](ParserBed::Data &i, const ParserProgress &)
    {
        cs.insert(i.cID);
        m2[i.cID].push_back(i);
    });

    FileName tmp = System::tmpFile();
    std::ofstream out(tmp);

    for (auto &c : cs)
    {
        for (const auto &l : Locus::inter<ParserBed::Data, Locus>(m1[c], m2[c]))
        {
            const auto i = std::find_if(m1[c].begin(), m1[c].end(), [&](const ParserBed::Data &x)
            {
                return x.l == l;
            });

            if (i != m1[c].end())
            {
                out << c << "\t" << l.start-1 << "\t" << l.end << "\t" << i->name << "\n";
            }
            else
            {
                const auto i = std::find_if(m1[c].begin(), m1[c].end(), [&](const ParserBed::Data &x)
                {
                    return x.l.overlap(l);
                });

                out << c << "\t" << l.start-1 << "\t" << l.end << "\t" << i->name << "\n";
            }
        }
    }
    
    out.close();
    return tmp;
}
