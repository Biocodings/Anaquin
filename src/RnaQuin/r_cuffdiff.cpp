/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "RnaQuin/r_cuffdiff.hpp"
#include "parsers/parser_cdiff.hpp"

using namespace Anaquin;

RCuffdiff::Stats RCuffdiff::stats(const FileName &file, const RCuffdiff::Options &o)
{
    Stats stats;
    
    ParserCDiff::parse(file, [&](const ParserCDiff::Data &x, const ParserProgress &)
    {
        if (x.status == ParserCDiff::Data::Status::Tested)
        {
            RCuffdiff::Stats::Data d;
            
            d.p    = x.p;
            d.q    = x.q;
            d.cID  = x.cID;
            d.gID  = x.gID;
            d.logF = x.logF;
            d.iID  = x.gID == x.tID ? "" : x.tID;

            stats.data.push_back(d);
        }
    });
    
    return stats;
}

void RCuffdiff::report(const FileName &file, const Options &o)
{
    const auto stats = RCuffdiff::stats(file, o);

    /*
     * Generating "RnaCuffdiff_converted.txt"
     */

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";

    o.generate("RnaCuffdiff_converted.txt");
    o.writer->open("RnaCuffdiff_converted.txt");

    o.writer->write((boost::format(format) % "ChrID"
                                           % "GeneID"
                                           % "IsoformID"
                                           % "FoldChange"
                                           % "FoldSE"
                                           % "PValue"
                                           % "QValue"
                                           % "Average").str());
    
    for (const auto &i : stats.data)
    {
        o.writer->write((boost::format(format) % i.cID
                                               % i.gID
                                               % (i.iID.empty() ? "-" : i.iID)
                                               % i.logF
                                               % "-"
                                               % p2str(i.p)
                                               % p2str(i.q)
                                               % "-").str());
    }
    
    o.writer->close();
}