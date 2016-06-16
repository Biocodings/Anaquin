#include "RnaQuin/r_cuffdiff.hpp"
#include "writers/file_writer.hpp"
#include "parsers/parser_cdiff.hpp"

using namespace Anaquin;

void RCuffdiff::analyze(const FileName &src, const FileName &output, const RCuffdiff::Options &o)
{
    FileWriter out(o.work);
    out.open(output);
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";
    out.write((boost::format(format) % "ChrID"
                                     % "GeneID"
                                     % "FoldChange"
                                     % "FoldSE"
                                     % "PValue"
                                     % "QValue"
                                     % "Average").str());
    
    ParserCDiff::parse(src, [&](const ParserCDiff::Data &x, const ParserProgress &)
    {
        out.write((boost::format(format) % x.cID
                                         % x.id
                                         % x.logF
                                         % "-"
                                         % x.p
                                         % x.q
                                         % "-").str());
    });
    
    out.close();
}

void RCuffdiff::report(const FileName &file, const Options &o)
{
    RCuffdiff::analyze(file, "RnaCuffdiff_converted.txt", o);
}