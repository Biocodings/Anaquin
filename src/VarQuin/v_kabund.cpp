#include "data/tokens.hpp"
#include "VarQuin/v_kabund.hpp"
#include "parsers/parser_salmon.hpp"

using namespace Anaquin;

extern Scripts PlotCNV();
extern Scripts PlotKAllele();
extern Scripts PlotConjoint();

typedef VKAbund::Software Software;

VKAbund::Stats VKAbund::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VKAbund::Stats stats;

    switch (o.soft)
    {
        case Software::Salmon:
        {
            ParserSalmon::parse(Reader(file), [&](const ParserSalmon::Data &x, const ParserProgress &)
            {
                auto found = false;
                
                switch (o.mode)
                {
                    case Mode::CNVLad:
                    {
                        auto n = noLast(x.name, "_");

                        if (r.match(n + "_CNV"))
                        {
                            stats.add(n + "_CNV", r.match(n + "_CNV")->concent(), x.abund);
                        }
                        
                        break;
                    }

                    case Mode::ConLad:
                    {
                        if (r.match(x.name + "_CON"))
                        {
                            stats.add(x.name + "_CON", r.match(x.name + "_CON")->concent(), x.abund);
                        }

                        break;
                    }
                        
                    case Mode::AFLad:
                    {
                        auto n = noLast(x.name, "_");
                     
                        if (r.match(n + "_AF"))
                        {
                            found = true;
                            stats.add(n + "_AF", r.match(n + "_AF")->concent(), x.abund);
                        }
                        
                        break;
                    }
                }

                if (!found)
                {
                    o.logInfo("Unknown: " + x.name);
                }
            });
        }
    }

    return stats;
}

static Scripts generateSummary(const FileName &src, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    extern FileName MixRef();
    
//    const auto &r = Standard::instance().r_var;
    const auto ls = stats.linear();
    
    const auto format = "-------VarKAbund Output Results\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------Detected sequins\n\n"
                        "       Synthetic: %4%\n"
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
                                  % "????" //r.countSeqs() // 2
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

static bool isCNV(const SequinID &x)    { return hasSub(x, "_CNV"); }
static bool isAllele(const SequinID &x) { return hasSub(x, "_AF");  }

static bool writeCNV(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    for (const auto &i : stats)
    {
        if (isCNV(i.first))
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
                    if (isStarted("GS", i.first) || isStarted("GI", i.first))
                    {
                        o.writer->write((boost::format(format) % noLast(i.first, "_")
                                                               % i.second.x
                                                               % i.second.y).str());
                    }
                }
            }
            
            o.writer->close();
            return true;
        }
    }
    
    return false;
}

static bool writeAllele(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
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
        
        if (isAllele(sID))
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
    
    if (!mr.empty())
    {
        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(format) % "Name"
                                               % "ExpRef"
                                               % "ExpVar"
                                               % "ObsRef"
                                               % "ObsVar").str());

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
    
    return !mr.empty();
}

static bool writeConjoint(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    for (const auto &i : stats)
    {
        if (!isAllele(i.first) && !isCNV(i.first))
        {
            o.generate(file);
            o.writer->open(file);
            o.writer->write((boost::format(format) % "Name" % "Sequin" % "Expected" % "Observed" % "Fold").str());
            
            for (const auto &i : stats)
            {
                if (!isAllele(i.first) && !isCNV(i.first))
                {
                    o.writer->write((boost::format(format) % i.first
                                                           % noLast(i.first, "_")
                                                           % i.second.x
                                                           % i.second.y
                                                           % last(i.first, "_")).str());
                }
            }
            
            o.writer->close();
            return true;
        }
    }
    
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

static void writeAlleleR(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
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

static void writeConjointR(const FileName &file, const VKAbund::Stats &stats, const VKAbund::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RWriter::createRLinear("VarKAbund_sequins.csv",
                                           o.work,
                                           "Expected Concentration vs Observed Abundance",
                                           "Expected Concentration (log2)",
                                           "Observed Abundance (log2)",
                                           "log2(data$Expected)",
                                           "log2(data$Observed)",
                                           "input",
                                           true,
                                           PlotConjoint()));
    o.writer->close();
}

void VKAbund::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);
    
    o.generate("VarKAbund_summary.stats");
    o.writer->open("VarKAbund_summary.stats");
    o.writer->write(generateSummary(file, stats, o));
    o.writer->close();

    switch (o.mode)
    {
        case Mode::CNVLad:
        {
            writeCNV("VarKAbund_sequins.csv", stats, o);
            writeCNVR("VarKAbund_linear.R",   stats, o);
            break;
        }

        case Mode::ConLad:
        {
            writeConjoint("VarKAbund_sequins.csv", stats, o);
            writeConjointR("VarKAbund_linear.R",   stats, o);
            break;
        }

        case Mode::AFLad:
        {
            writeAllele("VarKAbund_sequins.csv", stats, o);
            writeAlleleR("VarKAbund_linear.R",   stats, o);
            break;
        }
    }
}
