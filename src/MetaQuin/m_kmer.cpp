/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#include "MetaQuin/m_kmer.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMKMer();

// Defined by Kallisto
extern int __main__(int argc, char *argv[]);

MKMer::Stats MKMer::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MKMer::Stats stats;
    
    // Initialize the distribution for each sequin
    //stats.hist = r.hist();
    
    /*
     * Parsing the generated files. We're interested in the file listing the abundance.
     */
    
    ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &d, const ParserProgress &)
    {
        const auto m = r.match(d.id);
                              
        if (m)
        {
            // Expected abundance
            const auto known = m->mixes.at(Mix_1);
            
            // Measured abundance
            const auto measured = d.abund;
            
            if (measured)
            {
                stats.add(d.id, known, measured);
                
                stats.n_syn++;
                //stats.hist.at(d.id)++;
            }
        }
    });
    
    return stats;
}

void MKMer::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating MetaKmer_summary.stats
     */
    
    o.info("Generating MetaKmer_summary.stats");
    o.writer->open("MetaKmer_summary.stats");
//    o.writer->write(StatsWriter::inflectSummary(o.rAnnot,
//                                                o.rAnnot,
//                                                file,
//                                                stats.hist,
//                                                stats,
//                                                stats,
//                                                "sequins"));
    o.writer->close();
    
    /*
     * Generating MetaKmer_sequins.csv
     */
    
    o.info("Generating MetaKmer_sequins.csv");
    o.writer->open("MetaKmer_sequins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();
    
    /*
     * Generating MetaKmer_abund.R
     */
    
    o.info("Generating MetaKmer_abund.R");
    o.writer->open("MetaKmer_abund.R");
    o.writer->write(RWriter::createScript("MetaKmer_sequins.csv", PlotMKMer()));
    o.writer->close();
}