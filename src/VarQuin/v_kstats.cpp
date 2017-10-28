#include "Kallisto.hpp"
#include "ss/stats.hpp"
#include <boost/format.hpp>
#include "VarQuin/v_kstats.hpp"

using namespace Anaquin;

// Defined in Kallisto.cpp
extern KMStats Kallisto(const std::string &, const std::string &, const std::string &, unsigned);

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
    auto allIndex = o.allIndex;
    
    if (isEnd(allIndex, "fa") || isEnd(allIndex, "fasta"))
    {
//        std::cout << (boost::format("kallisto index -i /tmp/tmp.index %1%") % allIndex).str() << std::endl;
//        system((boost::format("kallisto index -i /tmp/tmp.index %1%") % allIndex).str().c_str());
//        allIndex = "/tmp/tmp.index";
    }

    // Run k-mer analysis
    const auto x = Kallisto(allIndex, files[0], files[1], 31);
    
    /*
     * 1: Dilution analysis
     */
    
    stats.nGen = x.nGen;
    stats.nSeq = x.nSeq;
    
    /*
     * 2: Allele frequency ladder (cancer sequins)
     */
    
    for (const auto &i : x.vars)
    {
        auto __count__ = [&](const std::vector<KMPair> &ps)
        {
            std::vector<unsigned> l;
            
            for (const auto &p : ps)
            {
                A_ASSERT(x.spans.count(p.normal));
                A_ASSERT(x.spans.count(p.revComp));
                
                const auto nc = x.spans.at(p.normal);
                const auto rc = x.spans.at(p.revComp);

                l.push_back(nc + rc);
            }
            
            return SS::median(l);
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

    o.writer->write((boost::format(format) % p1               // 1
                                           % p2               // 2
                                           % stats.nGen       // 3
                                           % stats.nSeq       // 4
                                           % stats.dilution() // 5
                                           % lm.m             // 6
                                           % lm.r             // 7
                                           % lm.R2            // 8
                                           % lm.F             // 9
                                           % lm.p             // 10
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

static void writeQuins(const FileName &file, const Stats &stats, const Options &o)
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
    
    o.info("Number of genome reads: " + std::to_string(stats.nGen));
    o.info("Number of sequin reads: " + std::to_string(stats.nSeq));
    
    /*
     * Generating VarKStats_summary.stats
     */
    
    writeSummary("VarKStats_summary.stats", files[0], files[1], stats, o);
    
    /*
     * Generating VarKStats_allele.tsv
     */
    
    writeQuins("VarKStats_allele.tsv", stats, o);
    
    /*
     * Generating VarKStats_allele.R
     */

    writeAlleleR("VarKStats_allele.R", "VarKStats_allele.tsv", stats, o);
}
