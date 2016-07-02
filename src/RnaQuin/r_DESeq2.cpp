#include "RnaQuin/r_DESeq2.hpp"
#include "writers/file_writer.hpp"
#include "parsers/parser_DESeq2.hpp"

using namespace Anaquin;

void RDESeq2::analyze(const FileName &src, const FileName &output, const RDESeq2::Options &o)
{
    FileWriter out(o.work);
    out.open(output);
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";
    out.write((boost::format(format) % "ChrID"
                                     % "GeneID"
                                     % "IsoformID"
                                     % "FoldChange"
                                     % "FoldSE"
                                     % "PValue"
                                     % "QValue"
                                     % "Average").str());

    ParserDESeq2::parse(src, [&](const ParserDESeq2::Data &x, const ParserProgress &)
    {
        out.write((boost::format(format) % x.cID
                                         % x.gID
                                         % "-"
                                         % x.logF
                                         % x.logFSE
                                         % x.p
                                         % x.q
                                         % x.baseMean).str());
    });

    out.close();
}

void RDESeq2::report(const FileName &file, const Options &o)
{
    RDESeq2::analyze(file, "RnaDESeq2_converted.txt", o);
}