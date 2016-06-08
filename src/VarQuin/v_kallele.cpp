#include "data/pachter.hpp"
#include "VarQuin/v_kallele.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotScatter();

VKAllele::Stats VKAllele::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    // Run quantification in Kallisto
    const auto abundFile = Pachter::externalQuant(o.index, file1, file2);

    VKAllele::Stats stats = VFreq::analyze(abundFile);
    
    return stats;
}

void VKAllele::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &stats = analyze(file1, file2, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating VarKAllele_summary.stats
     */

    o.info("Generating VarKAllele_summary.stats");
    o.writer->open("VarKAllele_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rAnnot,
                                                o.rAnnot,
                                                (file1 + " & " + file2),
                                                stats.hist,
                                                stats,
                                                stats.vars,
                                                "sequins"));
    o.writer->close();

    /*
     * Generating VarKAllele_quins.csv
     */

    o.info("Generating VarKAllele_quins.csv");
    o.writer->open("VarKAllele_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats.vars, "input", "measured"));
    o.writer->close();
    
    /*
     * Generating VarKAllele_allele.R
     */

    o.info("Generating VarKAllele_allele.R");
    o.writer->open("VarKAllele_allele.R");
    o.writer->write(RWriter::createScript("VarKAllele_quins.csv", PlotScatter()));
    o.writer->close();
}