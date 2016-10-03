#include <boost/format.hpp>
#include "data/standard.hpp"
#include "tools/gtf_data.hpp"
#include "RnaQuin/r_gene.hpp"

using namespace Anaquin;

RGene::Stats RGene::stats(const RGene::Options &o)
{
    const auto &r = Standard::instance().r_rna;
    
    RGene::Stats stats;
    
    for (const auto &gID : r.getGenes(ChrIS))
    {
        for (auto mix = 0; mix < r.countMix(); mix++)
        {
            auto &x = stats[static_cast<Mixture>(mix)][gID];
            
            x.gID = gID;
            x.len = r.findGene(ChrIS, gID)->l.length();
            x.con = r.concent(gID, static_cast<Mixture>(mix));
        }
    }
    
    return stats;
}

void RGene::report(const RGene::Options &o)
{
    const auto r = stats(o);
    o.writer->open("RnaGene_sequins.csv");

    if (r.size() == 1)
    {
        const auto format = "%1%\t%2%\t%3%";
        o.writer->write(((boost::format(format) % "ID" % "Length" % "Mix (attomol/ul)")).str());
        
        for (const auto &j : r.at(Mix_1))
        {
            o.writer->write(((boost::format(format) % j.second.gID % j.second.len % j.second.con)).str());
        }
    }
    else
    {
        const auto format = "%1%\t%2%\t%3%\t%4%";
        o.writer->write(((boost::format(format) % "ID" % "Length" % "Mix A (attomol/ul)" % "Mix B (attomol/ul)")).str());

        for (const auto &j : r.at(Mix_1))
        {
            o.writer->write(((boost::format(format) % j.second.gID
                                                    % j.second.len
                                                    % j.second.con
                                                    % r.at(Mix_2).at(j.second.gID).con)).str());
        }
    }
    
    o.writer->close();
}
