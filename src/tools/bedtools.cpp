#include <set>
#include <fstream>
#include "data/locus.hpp"
#include "tools/system.hpp"
#include "tools/bedtools.hpp"
#include "parsers/parser_bed.hpp"

using namespace Anaquin;

FileName BedTools::intersect(const FileName &x, const FileName &y)
{
    std::map<ChrID, std::vector<ParserBed::Data>> m;

    ParserBed::parse(x, [&](ParserBed::Data &i, const ParserProgress &) { m[i.cID].push_back(i); });
    ParserBed::parse(y, [&](ParserBed::Data &i, const ParserProgress &) { m[i.cID].push_back(i); });

    FileName tmp = System::tmpFile();
    std::ofstream out(tmp);

    for (auto &i : m)
    {
        for (const auto &l : Locus::merge<ParserBed::Data, Locus>(i.second))
        {
            out << i.first << "\t" << l.start-1 << "\t" << l.end << "\n";
        }
    }
    
    out.close();
    return tmp;
}
