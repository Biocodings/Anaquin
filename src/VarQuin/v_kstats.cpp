#include "Kallisto.hpp"
#include "ss/stats.hpp"
#include <boost/format.hpp>
#include "VarQuin/v_kstats.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotKAllele();

typedef VKStats::Stats Stats;
typedef VKStats::Options Options;

Stats VKStats::analyze(const std::vector<FileName> &files, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    // Ladder for allele frequency
    const auto l1 = r.seqsL1();

    Stats stats;
    o.logInfo("Generating a k-mer index for: " + o.fa);
    
    /*
     * 1: Generate a FASTA file for human regions
     */
    
    const auto hg = KHumanFASTA(o.fa);

    /*
     * 2: Generate temporary index
     *
     *     - Sequin sequences
     *     - Forward sequences
     *     - Spanning sequences
     */
    
    const auto i1 = KBuildIndex(o.fa, 31);
    const auto i2 = KBuildIndex(hg, 31);

    stats.kStats = KCount(i1, i2, files[0], files[1], o.k);

    /*
     * 3: Work out minimum, maximum and median statistics
     */
    
    auto &kStats = stats.kStats;
    
    // For all measured reference sequin k-mers...
    for (const auto &i : kStats.k2c)
    {
        // Make sure only either normal or reverse complement is counted but not both
        A_ASSERT(!kStats.k2c.count(revcomp(i.first)));
        
        // Eg: CI_019 has [10,1,5,6...]
        stats.s2c[KKM2Sequin(i.first, o.k)].push_back(i.second);
    }
    
    for (auto &i : stats.s2c)
    {
        stats.mins[i.first] = SS::min(i.second);
        stats.meds[i.first] = SS::med(i.second);
        stats.maxs[i.first] = SS::max(i.second);
    }
    
    A_ASSERT(!stats.mins.empty());
    A_ASSERT(!stats.meds.empty());
    A_ASSERT(!stats.maxs.empty());

    
//    for (const auto &i : kStats.vars)
//    {
//        auto __counts__ = [&](const SequinID &sID, const std::vector<KMPair> &x, std::vector<Counts> &v)
//        {
//            for (const auto &i : x)
//            {
//                const auto nc = kStats.k2c.count(i.norm) ? kStats.k2c.at(i.norm) : 0;
//                const auto rc = kStats.k2c.count(i.rcom) ? kStats.k2c.at(i.rcom) : 0;
//
//                // Don't bother if the k-mer is missing
//                if (nc + rc)
//                {
//                    v.push_back(nc + rc);
//                }
//            }
//        };
//
//        std::vector<Counts> v;
//        __counts__(i.first, i.second.R, v);
//        __counts__(i.first, i.second.V, v);
//
//        Counts min, med, max;
//
//        if (v.empty())
//        {
//            min = med = max = std::numeric_limits<Counts>::quiet_NaN();
//        }
//        else
//        {
//            min = SS::min(v);
//            med = SS::med(v);
//            max = SS::max(v);
//
//            A_ASSERT(min >= 0);
//            A_ASSERT(med >= 0);
//            A_ASSERT(max >= 0);
//        }
//
//        stats.mins[i.first] = min;
//        stats.meds[i.first] = med;
//        stats.maxs[i.first] = max;
//    }
    
    /*
     * 4: Allele frequency ladder
     */
    
    for (const auto &i : stats.kStats.vars)
    {
        auto __count__ = [&](const std::vector<KMPair> &ps)
        {
            std::vector<unsigned> l;
            
            for (const auto &p : ps)
            {
                A_ASSERT(stats.kStats.spans.count(p.norm));
                A_ASSERT(stats.kStats.spans.count(p.rcom));
                
                const auto nc = stats.kStats.spans.at(p.norm);
                const auto rc = stats.kStats.spans.at(p.rcom);

                l.push_back(nc + rc);
            }
            
            return SS::med(l);
        };
        
        const auto rn = __count__(i.second.R);
        const auto vn = __count__(i.second.V);
        
        if (rn + vn == 0)
        {
            o.logInfo(i.first + " not found");
            continue;
        }
        else if (vn == 0)
        {
            o.logInfo(i.first + " abundance is zero");
            continue;
        }
        
        // Expected allele frequency
        const auto exp = r.input1(i.first);
        
        // Measured allele frequency
        const auto obs = (float) vn / (rn + vn);
        
        // Allele frequency ladder
        stats.af.add(i.first, exp, obs);
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
                        "-------Number of reads\n\n"
                        "       Genome:   %3%\n"
                        "       Sequin:   %4%\n"
                        "       Dilution: %5%\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %6%\n"
                        "       Correlation: %7%\n"
                        "       R2:          %8%\n"
                        "       F-statistic: %9%\n"
                        "       P-value:     %10%\n";

    o.writer->write((boost::format(format) % p1                // 1
                                           % p2                // 2
                                           % stats.kStats.nGen // 3
                                           % stats.kStats.nSeq // 4
                                           % stats.dilution()  // 5
                                           % lm.m              // 6
                                           % lm.r              // 7
                                           % lm.R2             // 8
                                           % lm.F              // 9
                                           % lm.p              // 10
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
    const auto format = "%1%\t%2%\t%3%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Sequence"
                                           % "Counts").str());

    for (const auto &i : stats.kStats.k2c)
    {
        o.writer->write((boost::format(format) % KKM2Sequin(i.first, o.k)
                                               % i.first
                                               % i.second).str());
    }
        
    o.writer->close();
}

static void writeQuins(const FileName &file, const Stats &stats, const Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Minimum"
                                           % "Median"
                                           % "Maximum").str());
    
    for (const auto &i : stats.mins)
    {
        #define S(x) (std::isnan(x) ? "-" : std::to_string(x))
        
        o.writer->write((boost::format(format) % i.first
                                               % S(stats.mins.at(i.first))
                                               % S(stats.meds.at(i.first))
                                               % S(stats.maxs.at(i.first))).str());
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
    
    o.info("Number of genome reads: " + std::to_string(stats.kStats.nGen));
    o.info("Number of sequin reads: " + std::to_string(stats.kStats.nSeq));
    
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
