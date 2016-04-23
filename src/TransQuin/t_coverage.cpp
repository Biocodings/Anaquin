#include "TransQuin/t_coverage.hpp"

using namespace Anaquin;

TCoverage::Stats TCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance().r_trans;
    
    return CoverageTool::stats_(file, r.geneHist(ChrT), [&](const Alignment &align, const ParserProgress &)
    {
        if (align.cID == ChrT)
        {
            return r.findGene(ChrT, align.l, MatchRule::Contains);
        }

        return (const TransRef::GeneData *) nullptr;
    });
}

void TCoverage::report(const FileName &file, const TCoverage::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto stats = TCoverage::stats(file, o);

    /*
     * Generating summary statistics
     */
    
    o.info("Generating TransCoverage_summary.stats");
    
    CoverageTool::CoverageReportOptions x;
    
    x.rChrT   = o.rChrT;
    x.rGeno   = o.rGeno;
    x.writer  = o.writer;
    x.refs    = r.hist().size();
    x.length  = r.size();
    x.summary = "TransCoverage_summary.stats";
    
    CoverageTool::summary(stats, x, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        return r.match(Locus(i, j), MatchRule::Contains);
    });
    
    /*
     * Generating detailed CSV for the sequins
     */
    
    o.info("Generating VarCoverage_quins.csv");
    o.writer->open("VarCoverage_quins.csv");
    o.writer->write(CoverageTool::writeCSV(stats, x));
    o.writer->close();

    /*
     * Generating bedgraph for the standards
     */
    
    CoverageTool::CoverageBedGraphOptions y;
    
    y.writer = o.writer;
    y.file   = "TransCoverage_chrT.bedgraph";

    o.info("Generating TransCoverage_chrT.bedgraph");
    
    CoverageTool::bedGraph(stats, y, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        return r.findExon(ChrT, Locus(i, j), MatchRule::Contains);
    });
}