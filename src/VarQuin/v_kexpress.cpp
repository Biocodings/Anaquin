#include "data/tokens.hpp"
#include "VarQuin/v_kexpress.hpp"
#include "parsers/parser_salmon.hpp"

using namespace Anaquin;

extern Scripts PlotCNV();
extern Scripts PlotKLadder();
extern Scripts PlotKAllele();

typedef VKExpress::Software Software;

VKExpress::Stats VKExpress::analyze(const FileName &file, const FileName &cnv, const Options &o)
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
            
            ParserSalmon::parse(Reader(cnv), [&](const ParserSalmon::Data &x, const ParserProgress &)
            {
                auto a = noLast(x.name, "_");
                
                    if (r.match(a + "_CNV"))
                    {
                        stats.add(a + "_CNV", r.match(a + "_CNV")->concent(), x.abund);
                    }
            });
        }
    }

    return stats;
}

static Scripts generateSummary(const FileName &src, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
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
                                  % "????" //MixRef()      // 3
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

static bool isCancer(const std::string &x)
{
    return hasSub(x, "_R") || hasSub(x, "_V");
}

static bool isCNV(const SequinID &x) { return hasSub(x, "_CNV"); }

static void writeCNV(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    const auto format = "%1%\t%2%\t%3%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Expected"
                                           % "Observed").str());
    
    for (const auto &i : stats)
    {
        if (isCNV(i.first))
        {
            o.writer->write((boost::format(format) % noLast(i.first, "_")
                                                   % i.second.x
                                                   % i.second.y).str());
        }
    }
    
    o.writer->close();
}

static void writeCancer(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "ExpRef"
                                           % "ExpVar"
                                           % "ObsRef"
                                           % "ObsVar").str());
    
    std::set<SequinID> ids;
    std::map<SequinID, Concent> er, ev;
    std::map<SequinID, Measured> mr, mv;
    
    auto f = [&](const std::string &x)
    {
        return remove(remove(x, "_R"), "_V");
    };

    for (const auto &i : stats)
    {
        const auto &sID = i.first;
        
        if (isCancer(sID))
        {
            if (hasSub(sID, "_R"))
            {
                er[f(sID)] = i.second.x;
                mr[f(sID)] = i.second.y;
            }
            else
            {
                ev[f(sID)] = i.second.x;
                mv[f(sID)] = i.second.y;
            }
        }
    }

    for (const auto &i : mr)
    {
        const auto &id = i.first;

        o.writer->write((boost::format(format) % i.first
                                               % er.at(id)
                                               % ev.at(id)
                                               % (mr.count(id) ? mr.at(id) : 0)
                                               % (mv.count(id) ? mv.at(id) : 0)
                         ).str());
    }
    
    o.writer->close();
}

static void writeConjoint(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name" % "Sequin" % "Expected" % "Observed" % "Fold").str());

    for (const auto &i : stats)
    {
        if (!isCancer(i.first) && !isCNV(i.first))
        {
            o.writer->write((boost::format(format) % i.first
                                                   % noLast(i.first, "_")
                                                   % i.second.x
                                                   % i.second.y
                                                   % last(i.first, "_")).str());
        }
    }
    
    o.writer->close();
}

static void writeCNVR(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    extern std::string __full_command__;
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotCNV()) % date()
                                              % __full_command__
                                              % o.work
                                              % "VarKExpress_CNV.csv").str());
    o.writer->close();
}

static void writeCancerR(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    extern std::string __full_command__;

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotKAllele()) % date()
                                                  % __full_command__
                                                  % o.work
                                                  % "VarKExpress_cancer.csv").str());
    o.writer->close();
}

static void writeConjointR(const FileName &file, const VKExpress::Stats &stats, const VKExpress::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RWriter::createRLinear("VarKExpress_conjoint.csv",
                                           o.work,
                                           "Expected Concentration vs Observed Abundance",
                                           "Expected Concentration (log2)",
                                           "Observed Abundance (log2)",
                                           "log2(data$Expected)",
                                           "log2(data$Observed)",
                                           "input",
                                           true,
                                           PlotKLadder()));
    o.writer->close();
}

void VKExpress::report(const FileName &file, const FileName &file2, const Options &o)
{
    const auto stats = analyze(file, file2, o);
    
    /*
     * Generating VKExpress_summary.stats
     */
    
    o.generate("VarKExpress_summary.stats");
    o.writer->open("VarKExpress_summary.stats");
    o.writer->write(generateSummary(file, stats, o));
    o.writer->close();

    /*
     * Generating VarKExpress_CNV.csv
     */
    
    writeCNV("VarKExpress_CNV.csv", stats, o);
    
    /*
     * Generating VarKExpress_CNV.R
     */
    
    writeCNVR("VarKExpress_CNV.R", stats, o);
    
    /*
     * Generating VarKExpress_cancer.csv
     */
    
    writeCancer("VarKExpress_cancer.csv", stats, o);

    /*
     * Generating VarKExpress_conjoint.csv
     */
    
    writeConjoint("VarKExpress_conjoint.csv", stats, o);
    
    /*
     * Generating VarKExpress_cancer.R
     */
    
    writeCancerR("VarKExpress_cancer.R", stats, o);
    
    /*
     * Generating VarKExpress_conjoint.R
     */
    
    writeConjointR("VarKExpress_conjoint.R", stats, o);
}
