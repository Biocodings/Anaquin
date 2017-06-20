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

    // Ladder for conjoint
    const auto l2 = r.seqsL2();

    // Ladder for conjoint
    const auto l3 = r.seqsL3();

    VKAbund::Stats stats;

    ParserSalmon::parse(Reader(file), [&](const ParserSalmon::Data &x, const ParserProgress &)
    {
        /*
         * Ladder for allele frequency (homozygous and heterozygous for non-cancer)
         */
        
        if (l1.count(noLast(x.name, "_")))
        {
            if (x.name[x.name.size() - 1] == 'R')
            {
                stats.afR[noLast(x.name, "_")] = x.abund;
            }
            else
            {
                stats.afV[noLast(x.name, "_")] = x.abund;
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

static void writeQuins(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name" % "Copy" % "Allele" % "Ref" % "Var").str());

    for (const auto &i : stats.afV)
    {
        const auto &r = Standard::instance().r_var;
        o.writer->write((boost::format(format) % i.first
                                               % r.input4(i.first)
                                               % r.input1(i.first)
                                               % (stats.afR.count(i.first) ? stats.afR.at(i.first) : 0)
                                               % stats.afV.at(i.first)).str());
    }

    o.writer->close();
}

static bool writeConjoint(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
//    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
//    
//    for (const auto &i : stats)
//    {
//        if (!isAllele(i.first) && !isCNV(i.first))
//        {
//            o.generate(file);
//            o.writer->open(file);
//            o.writer->write((boost::format(format) % "Name" % "Sequin" % "Expected" % "Observed" % "Fold").str());
//            
//            for (const auto &i : stats)
//            {
//                if (!isAllele(i.first) && !isCNV(i.first))
//                {
//                    o.writer->write((boost::format(format) % i.first
//                                                           % noLast(i.first, "_")
//                                                           % i.second.x
//                                                           % i.second.y
//                                                           % last(i.first, "_")).str());
//                }
//            }
//            
//            o.writer->close();
//            return true;
//        }
//    }
    
    return false;
}

static void writeCNVR(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    extern std::string __full_command__;
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotCNV()) % date()
                                              % __full_command__
                                              % o.work
                                              % "VarKAbund_sequins.csv").str());
    o.writer->close();
}

static void writeAllele(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    extern std::string __full_command__;
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotKAllele()) % date()
                                                  % __full_command__
                                                  % o.work
                                                  % "VarKAbund_sequins.csv").str());
    o.writer->close();
}

//static void writeConjointR(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
//{
//    o.generate(file);
//    o.writer->open(file);
//    o.writer->write(RWriter::createRLinear("VarKAbund_sequins.csv",
//                                           o.work,
//                                           "Expected Concentration vs Observed Abundance",
//                                           "Expected Concentration (log2)",
//                                           "Observed Abundance (log2)",
//                                           "log2(data$Expected)",
//                                           "log2(data$Observed)",
//                                           "input",
//                                           true,
//                                           PlotConjoint()));
//    o.writer->close();
//}

void VKAbund::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarKAbund_summary.stats
     */
    
    writeSummary("VarKAbund_summary.stats", file, stats, o);

    /*
     * Generating VarKAbund_sequins.csv
     */
    
    writeQuins("VarKAbund_sequins.csv", stats, o);

    /*
     * Generating VarKAbund_allele.R
     */
    
    writeAllele("VarKAbund_allele.R", stats, o);
    
    /*
     * Generating VarKAbund_CNV.R
     */
    
    writeCNVR("VarKAbund_CNV.R", stats, o);

    
//            writeConjointR("VarKAbund_linear.R",   stats, o);
//            break;
}
