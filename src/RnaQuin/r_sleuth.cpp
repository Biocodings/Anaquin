#include "RnaQuin/r_sleuth.hpp"
#include "writers/file_writer.hpp"
#include "parsers/parser_sleuth.hpp"

using namespace Anaquin;

void RSleuth::analyze(const FileName &src, const FileName &output, const RSleuth::Options &o)
{
    FileWriter out(o.work);
    out.open(output);
    
    /*
     * Format: ChrID  Gene_ID  Fold_Change  Fold_SE  PValue  QValue  Average
     */

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";
    out.write((boost::format(format) % "ChrID"
                                     % "IsoID"
                                     % "FoldChange"
                                     % "FoldSE"
                                     % "PValue"
                                     % "QValue"
                                     % "Average").str());

    ParserSleuth::parse(src, [&](const ParserSleuth::Data &x, const ParserProgress &)
    {
        out.write((boost::format(format) % x.cID
                                         % x.id
                                         % x.logF
                                         % x.logFSE
                                         % x.p
                                         % x.q
                                         % x.baseMean).str());
    });

    out.close();
}

void RSleuth::report(const FileName &file, const Options &o)
{
    RSleuth::analyze(file, "RnaSleuth_converted.txt", o);
}