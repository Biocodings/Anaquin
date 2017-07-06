#include "data/tokens.hpp"
#include "VarQuin/v_kabund.hpp"
#include "parsers/parser_salmon.hpp"

using namespace Anaquin;

extern Scripts PlotCNV();
extern Scripts PlotKAllele();
extern Scripts PlotConjoint();

VKAbund::Stats VKAbund::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    // Ladder for allele frequency
    const auto l1 = r.seqsL1();

    // Ladder for conjoint (unit level)
    const auto l3 = r.seqsL3();

    // Ladder for germline mutations
    const auto l5 = r.seqsL5();
    
    A_ASSERT(!l3.empty());
    A_ASSERT(!l5.empty());
    
    VKAbund::Stats stats;

    ParserSalmon::parse(Reader(file), [&](const ParserSalmon::Data &x, const ParserProgress &)
    {
        /*
         * Ladder for germline mutation (0 and 0.5)
         */
        
        if (l5.count(noLast(x.name, "_")))
        {
            if (x.name[x.name.size() - 1] == 'R')
            {
                stats.germR[noLast(x.name, "_")] = x.abund;
            }
            else
            {
                stats.germV[noLast(x.name, "_")] = x.abund;
            }
        }
        
        /*
         * Ladder for conjoint sequins (unit level)
         */
        
        else if (l3.count(x.name))
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

static void writeSummary(const FileName &file, const FileName &src, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
//    o.generate("VarKAbund_summary.stats");
//    o.writer->open("VarKAbund_summary.stats");
//    o.writer->write(generateSummary(file, stats, o));
//    o.writer->close();
    
//    extern FileName MixRef();
//
////    const auto &r = Standard::instance().r_var;
//    const auto ls = stats.linear();
//
//    const auto format = "-------VarKAbund Output Results\n\n"
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

static void writeConjoint(const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate("VarKAbund_conjoint.csv");
    o.writer->open("VarKAbund_conjoint.csv");
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

    o.generate("VarKAbund_conjoint.R");
    o.writer->open("VarKAbund_conjoint.R");
    o.writer->write((boost::format(PlotConjoint()) % date()
                                                   % __full_command__
                                                   % o.work
                                                   % "VarKAbund_conjoint.csv"
                                                   % "Expected Copy Number vs Observed Abundance"
                                                   % "Expected Copy Number (log2)"
                                                   % "Observed Abundance (log2)"
                                                   % "log2(data$Input * data$Copy)"
                                                   % "log2(data$Observed)"
                     ).str());
    o.writer->close();
}

static void writeCNV(const VKAbund::Stats &stats, const VKAbund::Options &o)
{
//    extern std::string __full_command__;
//    
//    o.generate(file);
//    o.writer->open(file);
//    o.writer->write((boost::format(PlotCNV()) % date()
//                                              % __full_command__
//                                              % o.work
//                                              % "VarKAbund_sequins.csv").str());
//    o.writer->close();
}

static void writeSomatic(const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate("VarKAbund_somatic.csv");
    o.writer->open("VarKAbund_somatic.csv");
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
    
    o.generate("VarKAbund_somatic.R");
    o.writer->open("VarKAbund_somatic.R");
    o.writer->write((boost::format(PlotKAllele()) % date()
                                                  % __full_command__
                                                  % o.work
                                                  % "VarKAbund_somatic.csv").str());
    o.writer->close();
}

static void writeGerm(const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate("VarKAbund_germline.csv");
    o.writer->open("VarKAbund_germline.csv");
    o.writer->write((boost::format(format) % "Name"
                                           % "ExpFreq"
                                           % "ObsRef"
                                           % "ObsVar"
                                           % "ObsFreq").str());
    
    for (const auto &i : stats.germV)
    {
        const auto &r = Standard::instance().r_var;
        
        const auto obsR = (stats.germR.count(i.first) ? stats.germR.at(i.first) : 0);
        const auto obsV =  stats.germV.at(i.first);
        
        o.writer->write((boost::format(format) % i.first
                                               % r.input5(i.first)
                                               % obsR
                                               % obsV
                                               % ((float) obsV / (obsR + obsV))).str());
    }
    
    o.writer->close();
    
    extern std::string __full_command__;
    
    o.generate("VarKAbund_germline.R");
    o.writer->open("VarKAbund_germline.R");
    o.writer->write((boost::format(PlotKAllele()) % date()
                                                  % __full_command__
                                                  % o.work
                                                  % "VarKAbund_germline.csv").str());
    o.writer->close();
}

void VKAbund::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarKAbund_summary.stats
     */
    
    writeSummary("VarKAbund_summary.stats", file, stats, o);

    /*
     * Generating VarKAbund_somatic.csv and VarKAbund_somatic.R
     */
    
    writeSomatic(stats, o);

    /*
     * Generating VarKAbund_germline.csv and VarKAbund_germline.R
     */
    
    writeGerm(stats, o);
    
    /*
     * Generating VarKAbund_CNV.csv and VarKAbund_CNV.R
     */
    
    writeCNV(stats, o);

    /*
     * Generating VarKAbund_conjoint.csv and VarKAbund_conjoint.R
     */

    writeConjoint(stats, o);
}
