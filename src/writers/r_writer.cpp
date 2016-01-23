#include <iostream>
#include "stats/analyzer.hpp"
#include "data/experiment.hpp"
#include "writers/r_writer.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotLODR();

// Defined in resources.cpp
extern Scripts PlotMA();

// Defined in resources.cpp
extern Scripts PlotROC();

// Defined in resources.cpp
extern Scripts PlotSplice();

Scripts StatsWriter::inflectSummary()
{
    return "Summary for dataset: %1%\n\n"
           "   Synthetic:   %2% %3%\n\n"
           "   Experiment:  %4% %5%\n"
           "   Reference:   %6% %7%\n"
           "   Detected:    %8% %7%\n\n"
           "   ***\n"
           "   *** Detection Limits\n"
           "   ***\n\n"
           ""
           "   Break: %9% (%10%)\n\n"
           "   Left:  %11% + %12%x (R2 = %13%)\n"
           "   Right: %14% + %15%x (R2 = %16%)\n\n"
           "   ***\n"
           "   *** Statistics for linear regression\n"
           "   ***\n\n"
           "   Correlation: %17%\n"
           "   Slope:       %18%\n"
           "   R2:          %19%\n"
           "   F-statistic: %20%\n"
           "   P-value:     %21%\n"
           "   SSM:         %22%, DF: %23%\n"
           "   SSE:         %24%, DF: %25%\n"
           "   SST:         %26%, DF: %27%\n\n"
           "   ***\n"
           "   *** Statistics for linear regression (log2 scale)\n"
           "   ***\n\n"
           "   Correlation: %28%\n"
           "   Slope:       %29%\n"
           "   R2:          %30%\n"
           "   F-statistic: %31%\n"
           "   P-value:     %32%\n"
           "   SSM:         %33%, DF: %34%\n"
           "   SSE:         %35%, DF: %36%\n"
           "   SST:         %37%, DF: %38%\n";
}

Scripts StatsWriter::inflectSummary(const SInflectStats &stats)
{    
    return (boost::format(StatsWriter::inflectSummary()) % STRING(stats.files)      // 1
                                                         % STRING(stats.chrT_n)     // 2
                                                         % STRING(stats.chrT_p)     // 3
                                                         % STRING(stats.endo_n)     // 4
                                                         % STRING(stats.endo_p)     // 5
                                                         % STRING(stats.ref_n)      // 6
                                                         % STRING(stats.units)      // 7
                                                         % STRING(stats.det_n)      // 8
                                                         % STRING(stats.b)          // 9
                                                         % STRING(stats.bID)        // 10
                                                         % STRING(stats.lInt)       // 11
                                                         % STRING(stats.lSl)        // 12
                                                         % STRING(stats.lR2)        // 13
                                                         % STRING(stats.rInt)       // 14
                                                         % STRING(stats.rSl)        // 15
                                                         % STRING(stats.rR2)        // 16
                                                         % STRING(stats.nLog.r)     // 17
                                                         % STRING(stats.nLog.sl)    // 18
                                                         % STRING(stats.nLog.R2)    // 19
                                                         % STRING(stats.nLog.F)     // 20
                                                         % STRING(stats.nLog.p)     // 21
                                                         % STRING(stats.nLog.SSM)   // 22
                                                         % STRING(stats.nLog.SSM_D) // 23
                                                         % STRING(stats.nLog.SSE)   // 24
                                                         % STRING(stats.nLog.SSE_D) // 25
                                                         % STRING(stats.nLog.SST)   // 26
                                                         % STRING(stats.nLog.SST_D) // 27
                                                         % STRING(stats.wLog.r)
                                                         % STRING(stats.wLog.sl)
                                                         % STRING(stats.wLog.R2)
                                                         % STRING(stats.wLog.F)
                                                         % STRING(stats.wLog.p)
                                                         % STRING(stats.wLog.SSM)
                                                         % STRING(stats.wLog.SSM_D)
                                                         % STRING(stats.wLog.SSE)
                                                         % STRING(stats.wLog.SSE_D)
                                                         % STRING(stats.wLog.SST)
                                                         % STRING(stats.wLog.SST_D)
            ).str();
};

Scripts StatsWriter::inflectSummary(const std::vector<FileName> &files,
                                    const std::vector<LinearStats> &stats,
                                    const Units &units)
{
    SInflectStats r;
    r.units = units;
    
    for (auto i = 0; i < stats.size(); i++)
    {
        // Linear regression without logarithm
        const auto n_lm = stats[i].linear(false);
        
        // Linear regression with logarithm
        const auto l_lm = stats[i].linear(true);

        // Calcluate the inflect point after log-transformation
        const auto inf = stats[i].inflect(true);
        
        // Remember the break-point is on the log-scale, we'll need to convert it back
        const auto b = pow(2, inf.b);

        r.files.add(files[i]);

        r.b.add(b);
        r.bID.add(inf.id);
        r.lInt.add(inf.lInt);
        r.rInt.add(inf.rInt);
        r.lSl.add(inf.lSl);
        r.rSl.add(inf.rSl);
        
        r.lR2.add(inf.lR2);
        r.rR2.add(inf.rR2);

        r.nLog.p.add(n_lm.p);
        r.nLog.r.add(n_lm.r);
        r.nLog.F.add(n_lm.F);
        r.nLog.sl.add(n_lm.m);
        r.nLog.R2.add(n_lm.R2);
        r.nLog.SSM.add(n_lm.SSM);
        r.nLog.SSE.add(n_lm.SSE);
        r.nLog.SST.add(n_lm.SST);
        r.nLog.SSM_D.add(n_lm.SSM_D);
        r.nLog.SSE_D.add(n_lm.SSE_D);
        r.nLog.SST_D.add(n_lm.SST_D);

        r.wLog.p.add(l_lm.p);
        r.wLog.r.add(l_lm.r);
        r.wLog.F.add(l_lm.F);
        r.wLog.sl.add(l_lm.m);
        r.wLog.R2.add(l_lm.R2);
        r.wLog.SSM.add(l_lm.SSM);
        r.wLog.SSE.add(l_lm.SSE);
        r.wLog.SST.add(l_lm.SST);
        r.wLog.SSM_D.add(l_lm.SSM_D);
        r.wLog.SSE_D.add(l_lm.SSE_D);
        r.wLog.SST_D.add(l_lm.SST_D);
    }
    
    return inflectSummary(r);
}

Scripts RWriter::createSplice(const std::string &working, const FileName &fpkms)
{
    std::stringstream ss;
    ss << PlotSplice();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % working
                                    % fpkms).str();
}

Scripts RWriter::createLODR(const std::string &working, const FileName &dFile, const std::string &cFile, const std::string &lvl)
{
    std::stringstream ss;
    ss << PlotLODR();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % working
                                    % dFile
                                    % lvl).str();
}

Scripts RWriter::createMA(const std::string &working, const FileName &file, const std::string &lvl)
{
    std::stringstream ss;
    ss << PlotMA();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % working
                                    % file
                                    % lvl).str();
}

Scripts RWriter::createROC(const std::vector<FeatureID> &seqs, const std::vector<double> &ps, const std::string &lvl)
{
    assert(!seqs.empty() && seqs.size() == ps.size());

    std::stringstream ss;
    ss << PlotROC();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % ("\'" + boost::algorithm::join(seqs, "\',\'") + "\'")
                                    % concat(ps)
                                    % lvl).str();
}