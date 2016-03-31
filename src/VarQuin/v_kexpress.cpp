#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_kexpress.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVAbundAbund();

static void writeCSV(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    o.writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%";
    
    auto f = [&](const LinearStats &l)
    {
        const auto data = l.data(false);
        
        for (auto i = 0; i < data.ids.size(); i++)
        {
            o.writer->write((boost::format(format) % data.ids[i]
                                                   % data.x[i]
                                                   % data.y[i]).str());
        }
    };
    
    o.writer->write((boost::format(format) % "Sequin"
                                           % "EAbund"
                                           % "MAbund").str());
    f(stats);

    o.writer->close();
}

VKExpress::Stats VKExpress::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    
    
    
    
    
    VKExpress::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();

    return stats;
}

void VKExpress::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &stats = analyze(file1, file2, o);
    
    /*
     * Generating CSV for all sequins
     */
    
    writeCSV("VarKExpress_quins.csv", stats, o);
    
//    /*
//     * Generating for AlleleAllele
//     */
//    
//    o.writer->open("VarExpress_abundAbund.R");
//    o.writer->write(RWriter::createScript("VarExpress_quins.csv", PlotVAbundAbund()));
//    o.writer->close();
//    
//    /*
//     * Generating summary statistics
//     */
//    
//    o.writer->open("VarExpress_summary.stats");
//    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
//                                                o.rEndo,
//                                                file,
//                                                stats.hist,
//                                                stats,
//                                                stats,
//                                                "variants"));
//    o.writer->close();
}