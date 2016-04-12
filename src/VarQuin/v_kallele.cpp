#include "data/pachter.hpp"
#include "VarQuin/v_kallele.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVAllele();

VKAllele::Stats VKAllele::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VKAllele::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();
    
    stats.n_endo = NAN;
    
    // Run quantification in Kallisto
    const auto abundFile = Pachter::externalQuant(o.index, file1, file2);

    /*
     * Parsing the generated files. Obviosuly, we can't estimate the allele frequency unless we can detect both
     * reference and variant sequins.
     */

    std::set<SequinID> ids;
    std::map<SequinID, Coverage> matchr, matchv;

    ParserKallisto::parse(Reader(abundFile), [&](const ParserKallisto::Data &d, const ParserProgress &)
    {
        const auto match = r.match(d.id);
        
        if (match)
        {
            const auto bID = baseID(d.id);
            ids.insert(bID);
            
            if (isRefID(d.id))
            {
                matchr[bID] = d.abund;
            }
            else
            {
                matchv[bID] = d.abund;
            }
            
            stats.n_chrT++;
            stats.hist.at(d.id)++;
        }
    });

    for (const auto &id : ids)
    {
        if (matchr.count(id) && matchv.count(id))
        {
            const auto ref = matchr[id];
            const auto var = matchv[id];
            
            // Expected abundance
            const auto known = r.alleleFreq(id);
            
            // Measured abundance
            const auto measured = var / (ref + var);
            
            stats.add(id, known, measured);
        }
    }
    
    stats.limit = r.absolute(stats.hist);
    
    return stats;
}

void VKAllele::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &stats = analyze(file1, file2, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    o.info("Generating VarKAllele_summary.stats");
    o.writer->open("VarKAllele_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rEndo,
                                                (file1 + " & " + file2),
                                                stats.hist,
                                                stats,
                                                stats,
                                                "sequins"));
    o.writer->close();

    /*
     * Generating CSV for all sequins
     */

    o.info("Generating VarKAllele_quins.csv");
    o.writer->open("VarKAllele_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats, "EAlleleF", "MAlleleF"));
    o.writer->close();
    
    /*
     * Generating for allele vs allele
     */

    o.info("Generating VarKAllele_allele.R");
    o.writer->open("VarKAllele_allele.R");
    o.writer->write(RWriter::createScript("VarKAllele_quins.csv", PlotVAllele()));
    o.writer->close();
}