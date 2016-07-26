#include "RnaQuin/r_cufflink.hpp"
#include "writers/file_writer.hpp"
#include "parsers/parser_cufflink.hpp"

using namespace Anaquin;

void RCufflink::analyze(const FileName &src, const FileName &output, const RCufflink::Options &o)
{
    FileWriter out(o.work);
    out.open(output);
    
    const auto &r = Standard::instance().r_trans;

    Counts n_iso = 0;
    Counts countGen = 0;
    
    o.generate(output);
    
    ParserCufflink::parse(src, [&](const ParserCufflink::Data &x, const ParserProgress &)
    {
        if (Standard::isSynthetic(x.cID))
        {
            // Is this a reference isoform?
            if (r.match(x.tID))
            {
                n_iso++;
            }
            else if (r.findGene(x.cID, x.tID))
            {
                countGen++;
            }
        }
    });

    // Are we parsing genes.fpkm_tracking?
    const auto isGene = countGen > n_iso;
    
    /*
     * Format: ChrID  Gene_ID  Start  End  Abund
     */

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    if (isGene)
    {
        o.info("Gene Expression");
        out.write((boost::format(format) % "ChrID"
                                         % "GeneID"
                                         % "Start"
                                         % "End"
                                         % "Abund").str());
    }
    else
    {
        o.info("Isoform Expression");
        out.write((boost::format(format) % "ChrID"
                                         % "IsoID"
                                         % "Start"
                                         % "End"
                                         % "Abund").str());
    }

    ParserCufflink::parse(src, [&](const ParserCufflink::Data &x, const ParserProgress &)
    {
        if (isGene)
        {
            out.write((boost::format(format) % x.cID
                                             % x.tID
                                             % x.l.start
                                             % x.l.end
                                             % x.abund).str());
        }
        else
        {
            out.write((boost::format(format) % x.cID
                                             % x.tID
                                             % x.l.start
                                             % x.l.end
                                             % x.abund).str());
        }
    });

    out.close();
}

void RCufflink::report(const FileName &file, const Options &o)
{
    RCufflink::analyze(file, "RnaCufflink_converted.txt", o);
}