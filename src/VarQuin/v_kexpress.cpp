#include "tools/tools.hpp"
#include "VarQuin/v_kexpress.hpp"
#include "parsers/parser_salmon.hpp"

using namespace Anaquin;

extern Scripts PlotLinear_();

typedef VKExpress::Software Software;

VKExpress::Stats VKExpress::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VKExpress::Stats stats;

    switch (o.soft)
    {
        case Software::Salmon:
        {
            ParserSalmon::parse(Reader(file), [&](const ParserSalmon::Data &x, const ParserProgress &)
            {
                const auto m = r.match(x.name);
                
                if (m)
                {
                    stats.add(x.name, m->concent(), x.abund);
                }
            });
        }
    }

    return stats;
}

static Scripts generateSummary(const FileName &src, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    // Defined in resources.cpp
    extern FileName MixRef();
    
    const auto &r = Standard::instance().r_meta;
    const auto ls = stats.linear();
    
    const auto format = "-------VarKExpress Output\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference VarKExpress Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------Sequin Counts\n\n"
                        "       Synthetic: %4%\n"
                        "       Detection Sensitivity: %5% (attomol/ul) (%6%)\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %7%\n"
                        "       Correlation: %8%\n"
                        "       R2:          %9%\n"
                        "       F-statistic: %10%\n"
                        "       P-value:     %11%\n"
                        "       SSM:         %12%, DF: %13%\n"
                        "       SSE:         %14%, DF: %15%\n"
                        "       SST:         %16%, DF: %17%\n";
    
    const auto limit = stats.limitQuant();
    
    return (boost::format(format) % src           // 1
                                  % r.countSeqs() // 2
                                  % MixRef()      // 3
                                  % stats.size()  // 4
                                  % limit.abund   // 5
                                  % limit.id      // 6
                                  % ls.m          // 7
                                  % ls.r          // 8
                                  % ls.R2         // 9
                                  % ls.F          // 10
                                  % ls.p          // 11
                                  % ls.SSM        // 12
                                  % ls.SSM_D      // 13
                                  % ls.SSE        // 14
                                  % ls.SSE_D      // 15
                                  % ls.SST        // 16
                                  % ls.SST_D      // 17
            ).str();
}

static void writeQuins(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    const auto format = "%1%\t%2%\t%3%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name" % "Exp" % "Obs").str());

    for (const auto &i : stats)
    {
        o.writer->write((boost::format(format) % i.first % i.second.x % i.second.y).str());
    }
    
    o.writer->close();
}

static void writeRLinear(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RWriter::createRLinear("VKExpress_sequins.csv",
                                           o.work,
                                           "Expected Concentration vs Observed Abundance",
                                           "Expected Concentration (log2)",
                                           "Observed Abundance (log2)",
                                           "log2(data$Exp)",
                                           "log2(data$Obs)",
                                           "input",
                                           true,
                                           PlotLinear_()));
    o.writer->close();
}

void VKExpress::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
    /*
     * Generating VKExpress_summary.stats
     */
    
    o.generate("VKExpress_summary.stats");
    o.writer->open("VKExpress_summary.stats");
    o.writer->write(generateSummary(file, stats, o));
    o.writer->close();
    
    /*
     * Generating VKExpress_sequins.csv
     */
    
    writeQuins("VKExpress_sequins.csv", stats, o);
    
    /*
     * Generating VarKExpress_linear.R
     */
    
    writeRLinear("VarKExpress_linear.R", stats, o);
}
