#include "RnaQuin/r_cufflink.hpp"
#include "writers/file_writer.hpp"
#include "parsers/parser_cufflink.hpp"

using namespace Anaquin;

void RCufflink::analyze(const FileName &src, const FileName &output, const RCufflink::Options &o)
{
    FileWriter out(o.work);
    out.open(output);
    
    /*
     * Format: ChrID  Gene_ID  Start  End  Abund
     */

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    out.write((boost::format(format) % "ChrID"
                                     % "GeneID"
                                     % "Start"
                                     % "End"
                                     % "Abund").str());

    ParserCufflink::parse(src, [&](const ParserCufflink::Data &x, const ParserProgress &)
    {
        out.write((boost::format(format) % x.cID
                                         % x.id
                                         % x.l.start
                                         % x.l.end
                                         % x.abund).str());
    });

    out.close();
}

void RCufflink::report(const FileName &file, const Options &o)
{
    RCufflink::analyze(file, "RCufflink_converted.txt", o);
}