#include <map>
#include <iostream>
#include "standard.hpp"
#include "meta/m_blast.hpp"
#include "parsers/parser_blat.hpp"

using namespace Spike;

MBlast::Stats MBlast::analyze(const std::string &file)
{
    std::map<std::string, ParserBlat::BlatLine> psl;

    ParserBlat::parse(file, [&](const ParserBlat::BlatLine &l, const ParserProgress &)
    {
        psl[l.qName] = l;
    });

    /*
     * Create data-strucutre for each metaquin
     */
    
    std::map<MetaQuinID, MetaAlignment> m;

    const auto &mixB = Standard::instance().m_seq_B;
    
    for (const auto &seq : Standard::instance().m_seq_A)
    {
        m[seq.first].id   = seq.first;
        m[seq.first].mixA = seq.second.raw;
        m[seq.first].mixB = mixB.at(m[seq.first].id).abund();
    }
    
    /*
     * Compare each alignment to the metaquins
     */
    
    ParserBlat::parse(file, [&](const ParserBlat::BlatLine &l, const ParserProgress &)
    {
        if (m.count(l.tName))
        {
            m.at(l.tName).aligns.insert(Locus(l.tStart, l.tEnd));
        }
        else
        {
            std::cout << "Warning: " << l.tName << " not a metaquin" << std::endl;
        }
    });
    
    /*
     * Convert the results
     */
    
    MBlast::Stats stats;

    for (const auto &i : m)
    {
        stats.aligns.insert(i.second);
    }

    assert(stats.aligns.size() == m.size());
    
    return stats;
}