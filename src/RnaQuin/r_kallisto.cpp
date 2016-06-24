#include "data/standard.hpp"
#include "RnaQuin/r_kallisto.hpp"
#include "writers/file_writer.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

void RKallisto::analyze(const FileName &src, const FileName &output, const RKallisto::Options &o)
{
    FileWriter out(o.work);
    out.open(output);
    
    /*
     * Format: ChrID  Gene_ID  Start  End  Abund
     */
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    out.write((boost::format(format) % "ChrID"
                                     % "IsoID"
                                     % "Start"
                                     % "End"
                                     % "Abund").str());
    
    ParserKallisto::parse(src, [&](const ParserKallisto::Data &x, const ParserProgress &)
    {
        const auto &r = Standard::instance().r_trans;

        if (r.match(x.id))
        {
            out.write((boost::format(format) % "chrT"
                                             % x.id
                                             % "-"
                                             % "-"
                                             % x.abund).str());
        }
    });
    
    out.close();
}

void RKallisto::report(const FileName &file, const Options &o)
{
    RKallisto::analyze(file, "RnaKallisto_converted.txt", o);
}