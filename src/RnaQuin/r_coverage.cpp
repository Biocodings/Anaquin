#include "RnaQuin/r_coverage.hpp"

using namespace Anaquin;

RCoverage::Stats RCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance().r_trans;
    
    return CoverageTool::stats_(file, r.geneHist(ChrT), [&](const Alignment &align, const ParserProgress &)
    {
        if (align.cID == ChrT)
        {
            return r.findGene(ChrT, align.l, MatchRule::Contains);
        }

        throw "Not Implemented";
        //return (const TransRef::GeneData *) nullptr;
    });
}

void RCoverage::report(const FileName &file, const RCoverage::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto stats = RCoverage::stats(file, o);

    /*
     * Generating summary statistics
     */
    
    o.info("Generating RnaCoverage_summary.stats");
    
    CoverageTool::CoverageReportOptions x;
    
    x.rAnnot  = o.rAnnot;
    x.writer  = o.writer;
    x.refs    = r.hist().size();
    x.length  = r.size();
    x.summary = "RnaCoverage_summary.stats";
    
//    CoverageTool::summary(stats, x, [&](const ChrID &id, Base i, Base j, Coverage)
//    {
//        return r.match(Locus(i, j), MatchRule::Contains);
//    });
    
    /*
     * Generating RnaCoverage_sequins.csv
     */
    
    o.info("Generating RnaCoverage_sequins.csv");
    o.writer->open("RnaCoverage_sequins.csv");
    //o.writer->write(CoverageTool::writeCSV(stats, x));
    o.writer->close();

    /*
     * Generating RnaCoverage_sequins.bedgraph
     */
    
    CoverageTool::CoverageBedGraphOptions y;
    
    y.writer = o.writer;
    y.file   = "RnaCoverage_sequins.bedgraph";

    o.info("Generating RnaCoverage_sequins.bedgraph");
    
//    CoverageTool::bedGraph(stats, y, [&](const ChrID &id, Base i, Base j, Coverage)
//    {
//        return nullptr; // TODO: r.findExon(ChrT, Locus(i, j), MatchRule::Contains);
//    });
}