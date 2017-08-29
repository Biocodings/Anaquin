#include "VarQuin/v_kmer.hpp"
#include "parsers/parser_salmon.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

extern FileName AFRef();
extern Scripts  PlotKAllele();

static bool isKallisto(const FileName &file)
{
    return isEnded(file, "abundance.tsv");
}

VarKmer::Stats VarKmer::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    // Ladder for allele frequency
    const auto l1 = r.seqsL1();

    VarKmer::Stats stats;

    if (isKallisto(file))
    {
        ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &x, const ParserProgress &)
        {
            if (l1.count(noLast(x.id, "_")))
            {
                if (x.id[x.id.size() - 1] == 'R')
                {
                    stats.r[noLast(x.id, "_")] = x.abund;
                }
                else if (x.id[x.id.size() - 1] == 'V')
                {
                    stats.v[noLast(x.id, "_")] = x.abund;
                }
                else
                {
                    throw std::runtime_error("Unknown: " + x.id);
                }
            }
        });
    }
    
    for (const auto &i : stats.v)
    {
        if (stats.r.count(i.first))
        {
            stats.add(i.first, r.input1(i.first), ((float) i.second / (stats.r.at(i.first) + i.second)));
        }
        else
        {
            stats.add(i.first, r.input1(i.first), 1.0);
        }
    }
    
//    ParserSalmon::parse(Reader(file), [&](const ParserSalmon::Data &x, const ParserProgress &)
//    {
//        /*
//         * Ladder for conjoint sequins (unit level)
//         */
//        
//        if (l3.count(x.name))
//        {
//            stats.con[x.name] = x.abund;
//        }
//        
//        /*
//         * Ladder for allele frequency
//         */
//
//        else if (l1.count(noLast(x.name, "_")))
//        {
//            if (x.name[x.name.size() - 1] == 'R')
//            {
//                stats.canR[noLast(x.name, "_")] = x.abund;
//            }
//            else
//            {
//                stats.canV[noLast(x.name, "_")] = x.abund;
//            }
//        }
//    });

    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const VarKmer::Stats &stats, const VarKmer::Options &o)
{
    o.generate("VarKmer_summary.stats");
    o.writer->open("VarKmer_summary.stats");
    
    const auto &r = Standard::instance().r_var;
    
    // Log2 and skip zeros to keep results consistent with R
    const auto ls = stats.linear(true, true);

    const auto format = "-------VarKmer Output Results\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------Detected sequins\n\n"
                        "       Synthetic: %4%\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %7%\n"
                        "       Correlation: %8%\n"
                        "       R2:          %9%\n"
                        "       F-statistic: %10%\n"
                        "       P-value:     %11%\n";
    
    const auto limit = stats.limitQuant();

    o.writer->write((boost::format(format) % src               // 1
                                           % r.seqsL1().size() // 2
                                           % AFRef()           // 3
                                           % stats.size()      // 4
                                           % limit.abund       // 5
                                           % limit.id          // 6
                                           % ls.m              // 7
                                           % ls.r              // 8
                                           % ls.R2             // 9
                                           % ls.F              // 10
                                           % ls.p              // 11
                    ).str());
    o.writer->close();
}

static void writeQuins(const VarKmer::Stats &stats, const VarKmer::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate("VarKmer_quins.csv");
    o.writer->open("VarKmer_quins.csv");
    o.writer->write((boost::format(format) % "Name"
                                           % "ObsRef"
                                           % "ObsVar"
                                           % "ExpFreq"
                                           % "ObsFreq").str());
    
    for (const auto &i : stats)
    {
        const auto obsR = (stats.r.count(i.first) ? stats.r.at(i.first) : 0);
        const auto obsV =  stats.v.at(i.first);
        
        o.writer->write((boost::format(format) % i.first
                                               % obsR
                                               % obsV
                                               % i.second.x
                                               % i.second.y).str());
    }
    
    o.writer->close();

    extern std::string __full_command__;
    
    o.generate("VarKmer_linear.R");
    o.writer->open("VarKmer_linear.R");
    o.writer->write((boost::format(PlotKAllele()) % date()
                                                  % __full_command__
                                                  % o.work
                                                  % "VarKmer_quins.csv").str());
    o.writer->close();
}

void VarKmer::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarKmer_summary.stats
     */
    
    writeSummary("VarKmer_summary.stats", file, stats, o);

    /*
     * Generating VarKmer_quins.csv and VarKmer_linear.R
     */
    
    writeQuins(stats, o);
}
