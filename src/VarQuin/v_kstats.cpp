#include "Kallisto.hpp"
#include "ss/stats.hpp"
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_kstats.hpp"

using namespace Anaquin;

typedef VKStats::Stats Stats;
typedef VKStats::Options Options;

// Defined in resources.cpp
extern Scripts PlotKAllele();

Stats VKStats::analyze(const std::vector<FileName> &files, const Options &o)
{
    o.info("Threads: " + toString(o.thr));
    o.info("Index: " + o.fa);

    const auto &r = Standard::instance().r_var;
    
    // Ladder for allele frequency
    const auto l1 = r.seqsL1();

    Stats stats;
    
    o.logInfo("Generating index for: " + o.fa);
    const auto hg = KHumanFA(o.fa, stats.s2l);

    stats.kStats = KCount(KBuildIndex(o.fa, o.k),
                          KBuildIndex(hg, o.k),
                          files[0], files[1],
                          o.showReads ? o.work + "/VarKStats_genome.txt" : "",
                          o.showReads ? o.work + "/VarKStats_forward.txt" : "",
                          o.showReads ? o.work + "/VarKStats_sequins.txt" : "",
                          o.thr, o.k);
    auto &kStats = stats.kStats;

    /*
     * Work out k-mer statistics for each sequin
     */

    auto __stats__ = [&](const SequinCounts &x, Stats::Abundance &a)
    {
        for (const auto &i : x)
        {
            // Eg: CS_010
            const auto seq = noRV(i.first);
            
            for (const auto &j : i.second)
            {
                a.raws[seq].push_back(j.second);
            }
        }
    };
    
    __stats__(kStats.R.shared, stats.R);
    __stats__(kStats.R.uniqs,  stats.R);
    __stats__(kStats.F.shared, stats.F);
    __stats__(kStats.F.uniqs,  stats.F);

    if (stats.R.raws.empty() || stats.F.raws.empty())
    {
        o.warn("Sequin k-mers not found. Please check your input files.");
    }

    auto minMedMax = [&](Stats::Abundance &abund)
    {
        for (auto &i : abund.raws)
        {
            abund.sds[i.first]  = SS::SD(i.second);
            abund.mins[i.first] = SS::min(i.second);
            abund.meds[i.first] = SS::med(i.second);
            abund.maxs[i.first] = SS::max(i.second);
        }
    };
    
    minMedMax(stats.R);
    minMedMax(stats.F);

    if (stats.R.mins.empty() || stats.R.meds.empty() || stats.R.maxs.empty())
    {
        o.warn("Sequin k-mers not found. Please check your input files.");
    }
    
    /*
     * Allele frequency ladder
     */
    
    for (const auto &i : stats.kStats.seqs)
    {
        if (isCancer(i))
        {
            const auto R = i + "_R";
            const auto V = i + "_V";

            const auto rn = kStats.R.uniqs.count(R) ? SS::med_(kStats.R.uniqs.at(R)) : 0;
            const auto vn = kStats.R.uniqs.count(V) ? SS::med_(kStats.R.uniqs.at(V)) : 0;

            if (rn + vn == 0)
            {
                o.logInfo(R + " not found");
                continue;
            }
            else if (vn == 0)
            {
                o.logInfo(V + " abundance is zero");
                continue;
            }

            // Expected allele frequency
            const auto exp = r.input1(i);
            
            // Measured allele frequency
            const auto obs = (float) vn / (rn + vn);
            
            // Allele frequency ladder
            stats.af.add(i, exp, obs);
        }
    }

    return stats;
}

static void writeSummary(const FileName &file, const FileName &p1, const FileName &p2, const Stats &stats, const VKStats::Options &o)
{
    o.generate(file);
    o.writer->open(file);

    LinearModel lm;
    
    try { lm = stats.af.linear(true, true); }
    catch(...) {}

    const auto format = "-------VarKStats Output Results\n\n"
                        "       Summary for input: %1% and %2%\n\n"
                        "-------Alignment reads\n\n"
                        "       Genome:      %3%\n"
                        "       Sequin:      %4%\n"
                        "       Dilution:    %5%\n"
                        "       Error Rate:  %6%\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %7%\n"
                        "       Correlation: %8%\n"
                        "       R2:          %9%\n"
                        "       F-statistic: %10%\n"
                        "       P-value:     %11%\n";

    o.writer->write((boost::format(format) % p1                     // 1
                                           % p2                     // 2
                                           % stats.kStats.R.nMatch  // 3
                                           % stats.kStats.R.nNMatch // 4
                                           % stats.dilution()       // 5
                                           % stats.error()          // 6
                                           % lm.m                   // 7
                                           % lm.r                   // 8
                                           % lm.R2                  // 9
                                           % lm.F                   // 10
                                           % lm.p                   // 11
                    ).str());
    o.writer->close();
}

static void writeAlleleR(const FileName &file, const FileName &src, const Stats &stats, const Options &o)
{
    extern std::string __full_command__;
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotKAllele()) % date()
                                                  % __full_command__
                                                  % o.work
                                                  % src).str());
    o.writer->close();
}

static void writeKmers(const FileName &file, const Stats &stats, const Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Sequence"
                                           % "Counts"
                                           % "Type").str());

    for (const auto &i : stats.kStats.R.shared)
    {
        for (const auto &j : i.second)
        {
            o.writer->write((boost::format(format) % i.first
                                                   % j.first
                                                   % j.second
                                                   % "Shared").str());
        }
    }

    for (const auto &i : stats.kStats.R.uniqs)
    {
        for (const auto &j : i.second)
        {
            o.writer->write((boost::format(format) % i.first
                                                   % j.first
                                                   % j.second
                                                   % "Unique").str());
        }
    }

    o.writer->close();
}

static void writeQuins(const FileName &file, const Stats &stats, const Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Minimum (Reverse)"
                                           % "Median (Reverse)"
                                           % "Maximum (Reverse)"
                                           % "SD (Reverse)"
                                           % "Minimum (Forward)"
                                           % "Median (Forward)"
                                           % "Maximum (Forward)"
                                           % "SD (Forward)").str());
    
    for (const auto &seq : stats.kStats.seqs)
    {
        #define S(x,y) (x.count(y) ? std::isnan(x.at(y)) ? "-" : std::to_string(x.at(y)) : "-")
        
        o.writer->write((boost::format(format) % seq
                                               % S(stats.R.mins, seq)
                                               % S(stats.R.meds, seq)
                                               % S(stats.R.maxs, seq)
                                               % S(stats.R.sds,  seq)
                                               % S(stats.F.mins, seq)
                                               % S(stats.F.meds, seq)
                                               % S(stats.F.maxs, seq)
                                               % S(stats.F.sds,  seq)).str());
    }
    
    o.writer->close();
}

static void writeAllele(const FileName &file, const Stats &stats, const Options &o)
{
    const auto format = "%1%\t%2%\t%3%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "ExpFreq"
                                           % "ObsFreq").str());
    
    for (const auto &i : stats.af)
    {
        o.writer->write((boost::format(format) % i.first
                                               % i.second.x
                                               % i.second.y).str());
    }
    
    o.writer->close();
}

void VKStats::report(const std::vector<FileName> &files, const Options &o)
{
    A_ASSERT(files.size() == 2);    
    const auto stats = analyze(files, o);
    
    /*
     * Generating VarKStats_summary.stats
     */
    
    writeSummary("VarKStats_summary.stats", files[0], files[1], stats, o);
    
    /*
     * Genetating VarKStats_kmers.tsv
     */

    writeKmers("VarKStats_kmers.tsv", stats, o);
    
    /*
     * Genetating VarKStats_sequins.tsv
     */

    writeQuins("VarKStats_sequins.tsv", stats, o);

    /*
     * Generating VarKStats_allele.tsv
     */
    
    writeAllele("VarKStats_allele.tsv", stats, o);
    
    /*
     * Generating VarKStats_allele.R
     */

    writeAlleleR("VarKStats_allele.R", "VarKStats_allele.tsv", stats, o);
}
