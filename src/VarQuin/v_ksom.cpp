#include "data/tokens.hpp"
#include "VarQuin/v_ksom.hpp"
#include "parsers/parser_salmon.hpp"

using namespace Anaquin;

extern Scripts PlotCNV();
extern Scripts PlotKAllele();
extern Scripts PlotConjoint();

VarKSomatic::Stats VarKSomatic::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    // Ladder for allele frequency
    const auto l1 = r.seqsL1();

    // Ladder for conjoint (unit level)
    const auto l3 = r.seqsL3();

    A_ASSERT(!l3.empty());
    
    VarKSomatic::Stats stats;

    ParserSalmon::parse(Reader(file), [&](const ParserSalmon::Data &x, const ParserProgress &)
    {
        /*
         * Ladder for conjoint sequins (unit level)
         */
        
        if (l3.count(x.name))
        {
            stats.con[x.name] = x.abund;
        }
        
        /*
         * Ladder for allele frequency
         */

        else if (l1.count(noLast(x.name, "_")))
        {
            if (x.name[x.name.size() - 1] == 'R')
            {
                stats.canR[noLast(x.name, "_")] = x.abund;
            }
            else
            {
                stats.canV[noLast(x.name, "_")] = x.abund;
            }
        }
    });

    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const VarKSomatic::Stats &stats, const VarKSomatic::Options &o)
{
//    o.generate("VarKSomatic_summary.stats");
//    o.writer->open("VarKSomatic_summary.stats");
//    o.writer->write(generateSummary(file, stats, o));
//    o.writer->close();
    
//    extern FileName MixRef();
//
////    const auto &r = Standard::instance().r_var;
//    const auto ls = stats.linear();
//
//    const auto format = "-------VarKSomatic Output Results\n\n"
//                        "       Summary for input: %1%\n\n"
//                        "-------Reference Annotations\n\n"
//                        "       Synthetic: %2%\n"
//                        "       Mixture file: %3%\n\n"
//                        "-------Detected sequins\n\n"
//                        "       Synthetic: %4%\n"
//                        "-------Linear regression (log2 scale)\n\n"
//                        "       Slope:       %7%\n"
//                        "       Correlation: %8%\n"
//                        "       R2:          %9%\n"
//                        "       F-statistic: %10%\n"
//                        "       P-value:     %11%\n"
//                        "       SSM:         %12%, DF: %13%\n"
//                        "       SSE:         %14%, DF: %15%\n"
//                        "       SST:         %16%, DF: %17%\n";
//    
//    const auto limit = stats.limitQuant();
//    
//    return (boost::format(format) % src           // 1
//                                  % "????" //r.countSeqs() // 2
//                                  % MixRef()      // 3
//                                  % stats.size()  // 4
//                                  % limit.abund   // 5
//                                  % limit.id      // 6
//                                  % ls.m          // 7
//                                  % ls.r          // 8
//                                  % ls.R2         // 9
//                                  % ls.F          // 10
//                                  % ls.p          // 11
//                                  % ls.SSM        // 12
//                                  % ls.SSM_D      // 13
//                                  % ls.SSE        // 14
//                                  % ls.SSE_D      // 15
//                                  % ls.SST        // 16
//                                  % ls.SST_D      // 17
//            ).str();
}

static void writeConjoint(const VarKSomatic::Stats &stats, const VarKSomatic::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate("VarKSomatic_conjoint.csv");
    o.writer->open("VarKSomatic_conjoint.csv");
    o.writer->write((boost::format(format) % "Name"
                                           % "Unit"
                                           % "Input"
                                           % "Copy"
                                           % "Observed").str());

    for (const auto &i : stats.con)
    {
        const auto &r = Standard::instance().r_var;
        
        o.writer->write((boost::format(format) % noLast(i.first, "_")
                                               % i.first
                                               % r.input2(noLast(i.first, "_"))
                                               % r.input3(i.first)
                                               % i.second).str());
    }
    
    o.writer->close();

    o.generate("VarKSomatic_conjoint.R");
    o.writer->open("VarKSomatic_conjoint.R");
    o.writer->write((boost::format(PlotConjoint()) % date()
                                                   % __full_command__
                                                   % o.work
                                                   % "VarKSomatic_conjoint.csv"
                                                   % "Expected Copy Number vs Observed Abundance"
                                                   % "Expected Copy Number (log2)"
                                                   % "Observed Abundance (log2)"
                                                   % "log2(data$Input * data$Copy)"
                                                   % "log2(data$Observed)"
                     ).str());
    o.writer->close();
}

static void writeCNV(const VarKSomatic::Stats &stats, const VarKSomatic::Options &o)
{
//    extern std::string __full_command__;
//    
//    o.generate(file);
//    o.writer->open(file);
//    o.writer->write((boost::format(PlotCNV()) % date()
//                                              % __full_command__
//                                              % o.work
//                                              % "VarKSomatic_sequins.csv").str());
//    o.writer->close();
}

static void writeSomatic(const VarKSomatic::Stats &stats, const VarKSomatic::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate("VarKSomatic_somatic.csv");
    o.writer->open("VarKSomatic_somatic.csv");
    o.writer->write((boost::format(format) % "Name"
                                           % "ExpFreq"
                                           % "ObsRef"
                                           % "ObsVar"
                                           % "ObsFreq").str());
    
    for (const auto &i : stats.canV)
    {
        const auto &r = Standard::instance().r_var;
        
        const auto obsR = (stats.canR.count(i.first) ? stats.canR.at(i.first) : 0);
        const auto obsV =  stats.canV.at(i.first);
        
        o.writer->write((boost::format(format) % i.first
                                               % r.input1(i.first)
                                               % obsR
                                               % obsV
                                               % ((float) obsV / (obsR + obsV))).str());
    }
    
    o.writer->close();

    extern std::string __full_command__;
    
    o.generate("VarKSomatic_somatic.R");
    o.writer->open("VarKSomatic_somatic.R");
    o.writer->write((boost::format(PlotKAllele()) % date()
                                                  % __full_command__
                                                  % o.work
                                                  % "VarKSomatic_somatic.csv").str());
    o.writer->close();
}

void VarKSomatic::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarKSomatic_summary.stats
     */
    
    writeSummary("VarKSomatic_summary.stats", file, stats, o);

    /*
     * Generating VarKSomatic_somatic.csv and VarKSomatic_somatic.R
     */
    
    writeSomatic(stats, o);

    /*
     * Generating VarKSomatic_CNV.csv and VarKSomatic_CNV.R
     */
    
    writeCNV(stats, o);

    /*
     * Generating VarKSomatic_conjoint.csv and VarKSomatic_conjoint.R
     */

    writeConjoint(stats, o);
}
