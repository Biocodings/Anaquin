#include <thread>
#include "data/standard.hpp"
#include "MetaQuin/m_diff.hpp"

using namespace Anaquin;

MDiff::Stats MDiff::analyze(const std::vector<FileName> &files, const Options &o)
{
    const auto &r = Standard::instance().r_meta;

    MDiff::Stats stats;

    MAbund::Options ao;
    ao.format = o.format;
    ao.writer = o.writer;

    std::thread t1([&]()
    {
        switch (o.format)
        {
            case Format::BAM:
            {
                stats.stats1 = MAbund::analyze(std::vector<FileName> { files[0] }, ao);
                break;
            }

            case Format::RayMeta:
            {
                stats.stats1 = MAbund::analyze(std::vector<FileName> { files[0], files[1] }, ao);
                break;
            }
        }
    });

    std::thread t2([&]()
    {
        switch (o.format)
        {
            case Format::BAM:
            {
                stats.stats1 = MAbund::analyze(std::vector<FileName> { files[1] }, ao);
                break;
            }

            case Format::RayMeta:
            {
                stats.stats2 = MAbund::analyze(std::vector<FileName> { files[2], files[3] }, ao);
                break;
            }
        }
    });
    
    t1.join();
    t2.join();

    std::map<SequinID, Measured> s1, s2;
    
    for (const auto &i : stats.stats1)
    {
        s1[i.first] = i.second.y;
    }

    for (const auto &i : stats.stats2)
    {
        s2[i.first] = i.second.y;
    }
    
    for (const auto &i : s1)
    {
        if (s2.count(i.first))
        {
            if (r.match(i.first))
            {
                stats.add(i.first, r.match(i.first)->fold(), s2[i.first] / i.second);
            }
            else
            {
                o.logWarn(i.first + " is not defined in the reference mixture file");
            }
        }
    }
    
    return stats;
}

static void writeQuins(const FileName &file, const MDiff::Stats &stats, const MDiff::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    const auto format = "%1%\t%2%\t%3%\t%4%";

    o.generate(file);
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID" % "Length" % "ExpFold" % "ObsFold").str());
    
    for (const auto &i : stats)
    {
        o.writer->write((boost::format(format) % i.first
                                               % r.match(i.first)->l.length()
                                               % i.second.x
                                               % i.second.y).str());
    }
    
    o.writer->close();
}

Scripts MDiff::generateRLinear(const FileName &src, const Stats &stats, const Options &o)
{
    return RWriter::createRLinear(src,
                                  o.work,
                                  "Differential Fold",
                                  "Expected Fold (log2)",
                                  "Measured Fold (log2)",
                                  "log2(data$ExpFold)",
                                  "log2(data$ObsFold)",
                                  "input",
                                  true);
}

static void writeRLinear(const FileName &file, const FileName &src, const MDiff::Stats &stats, const MDiff::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(MDiff::generateRLinear(src, stats, o));
    o.writer->close();
}

static Scripts generateSummary(const FileName &f1, const FileName &f2, const MDiff::Stats &stats, const MDiff::Options &o)
{
    // Defined in resources.cpp
    extern FileName MixRef();
    
    const auto &r = Standard::instance().r_meta;
    const auto ls = stats.linear();
    
    const auto format = "-------MetaFoldChange Output\n\n"
                        "       Summary for input: %1% and %2%\n\n"
                        "-------Reference MetaQuin Annotations\n\n"
                        "       Synthetic: %3%\n"
                        "       Mixture file: %4%\n\n"
                        "-------Sequin Counts\n\n"
                        "       Synthetic: %5%\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %6%\n"
                        "       Correlation: %7%\n"
                        "       R2:          %8%\n"
                        "       F-statistic: %9%\n"
                        "       P-value:     %10%\n"
                        "       SSM:         %11%, DF: %12%\n"
                        "       SSE:         %13%, DF: %14%\n"
                        "       SST:         %15%, DF: %16%\n";
    
    return (boost::format(format) % f1            // 1
                                  % f2            // 2
                                  % "????" //r.countSeqs() // 3
                                  % MixRef()      // 4
                                  % stats.size()  // 5
                                  % ls.m          // 6
                                  % ls.r          // 7
                                  % ls.R2         // 8
                                  % ls.F          // 9
                                  % ls.p          // 10
                                  % ls.SSM        // 11
                                  % ls.SSM_D      // 12
                                  % ls.SSE        // 13
                                  % ls.SSE_D      // 14
                                  % ls.SST        // 15
                                  % ls.SST_D      // 16
            ).str();
}

void MDiff::report(const std::vector<FileName> &files, const MDiff::Options &o)
{
    const auto stats = MDiff::analyze(files, o);
    
    /*
     * Generating MetaFoldChange_summary.stats
     */
    
    o.generate("MetaFoldChange_summary.stats");
    o.writer->open("MetaFoldChange_summary.stats");
    
    switch (o.format)
    {
        case Format::BAM:
        {
            o.writer->write(generateSummary(files[0], files[1], stats, o));
            break;
        }

        case Format::RayMeta:
        {
            o.writer->write(generateSummary(files[0], files[2], stats, o));
            break;
        }
    }

    o.writer->close();
    
    /*
     * Generating MetaFoldChange_sequins.csv
     */

    writeQuins("MetaFoldChange_sequins.csv", stats, o);
    
    /*
     * Generating MetaFoldChange_fold.R
     */
    
    writeRLinear("MetaFoldChange_fold.R", "MetaFoldChange_sequins.csv", stats, o);
}
