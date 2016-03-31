#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_kexpress.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVAbundAbund();

// Defined by Kallisto
extern int __main__(int argc, char *argv[]);

VKExpress::Stats VKExpress::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VKExpress::Stats stats;
    
    char *argv[8];
    
    const auto temp = std::string("/tmp/output");
    
    /*
     * Eg: kallisto quant -i /data/index/VarQuin/AVA010.v032.index -o output LVA086.1_val_1.fq LVA086.2_val_2.fq
     */
    
    argv[0] = new char[9];
    argv[1] = new char[6];
    argv[2] = new char[3];
    argv[3] = new char[o.file.size()+1];
    argv[4] = new char[3];
    argv[5] = new char[temp.size()+1];
    argv[6] = new char[file1.size()+1];
    argv[7] = new char[file2.size()+1];

    strcpy(argv[0], "kallisto");
    strcpy(argv[1], "quant");
    strcpy(argv[2], "-i");
    strcpy(argv[3], o.file.c_str());
    strcpy(argv[4], "-o");
    strcpy(argv[5], temp.c_str());
    strcpy(argv[6], file1.c_str());
    strcpy(argv[7], file2.c_str());

    optind = optreset = 1;
    
    /*
     * Execute Kallisto as if it were a standalone program
     */

    //__main__(8, argv);
    
    /*
     * Parse the generated files. We're interested in the file listing the abundance.
     */
    
    ParserKallisto::parse(Reader(temp + "/abundance.tsv"), [&](const ParserKallisto::Data &d, const ParserProgress &)
    {
        const auto m = r.match(d.id);
        
        if (m)
        {
            // Expected abundance
            const auto known = m->mixes.at(Mix_1);
            
            // Measured abundance
            const auto measured = d.abund;
            
            stats.add(d.id, known, measured);
        }
    });
    
    return stats;
}

void VKExpress::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &stats = analyze(file1, file2, o);

    /*
     * Generating summary statistics
     */

    o.writer->open("VarKExpress_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rEndo,
                                                (file1 + " & " + file2),
                                                stats.hist,
                                                stats,
                                                stats,
                                                "variants"));
    o.writer->close();
    
    /*
     * Generating CSV for all sequins
     */
    
    o.writer->open("VarKExpress_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating for AbundAbund
     */
    
    o.writer->open("VarKExpress_abundAbund.R");
    o.writer->write(RWriter::createScript("VarKExpress_quins.csv", PlotVAbundAbund()));
    o.writer->close();
}