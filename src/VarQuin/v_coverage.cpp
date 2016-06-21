#include "data/path.hpp"
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_coverage.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVDensity();

// Defined in resources.cpp
extern FileName refFile();

VCoverage::Stats VCoverage::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    throw "Not Implemented";
    
//    const auto &r = Standard::instance().r_var;
//
//    return CoverageTool::stats_(file, r.baseHist(), [&](const Alignment &align, const ParserProgress &)
//    {
//        if (align.cID == ChrT)
//        {
//            const auto m = r.match(align.l, MatchRule::Contains);
//            
//            if (m)
//            {
//                return r.findGene(baseID(m->id));
//            }
//        }
//
//        return (const VarRef::Base *) nullptr;
//    });
}

void VCoverage::report(const FileName &file, const VCoverage::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto stats = VCoverage::analyze(file, o);

    /*
     * Generating VarCoverage_summary.stats
     */
    
    o.generate("VarCoverage_summary.stats");

    CoverageTool::CoverageReportOptions x;
    
    x.rAnnot   = o.rAnnot;
    //x.rGeno   = o.rGeno;
    x.writer  = o.writer;
    x.refs    = r.hist().size();
    x.length  = r.size();
    x.summary = "VarCoverage_summary.stats";
    
//    CoverageTool::summary(stats, x, [&](const ChrID &id, Base i, Base j, Coverage)
//    {
//        return r.match(Locus(i, j), MatchRule::Contains);
//    });

    /*
     * Generating VarCoverage_sequins.csv
     */

    o.generate("VarCoverage_sequins.csv");
    o.writer->open("VarCoverage_sequins.csv");
    //o.writer->write(CoverageTool::writeCSV(stats, x));
    o.writer->close();

    /*
     * Generating VarCoverage_quins.bedgraph
     */

    CoverageTool::CoverageBedGraphOptions y;
    
    y.writer = o.writer;
    y.file   = "VarCoverage_quins.bedgraph";

    o.generate("VarCoverage_quins.bedgraph");

//    CoverageTool::bedGraph(stats, y, [&](const ChrID &id, Base i, Base j, Coverage)
//    {
//        return r.match(Locus(i, j), MatchRule::Contains);
//    });

    /*
     * Generating VarCoverage_density.R
     */

    o.generate("VarCoverage_density.R");
    o.writer->open("VarCoverage_density.R");
    o.writer->write(RWriter::createScript("VarCoverage_quins.bedgraph", PlotVDensity(), path2file(refFile())));
    o.writer->close();
}