#include <iostream>
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
                                                         % STRING(stats.lInter)     // 11
                                                         % STRING(stats.lSlope)     // 12
                                                         % STRING(stats.lR2)        // 13
                                                         % STRING(stats.rInter)     // 14
                                                         % STRING(stats.rSlope)     // 15
                                                         % STRING(stats.rR2)        // 16
                                                         % STRING(stats.nLog.cor)   // 17
                                                         % STRING(stats.nLog.slope) // 18
                                                         % STRING(stats.nLog.R2)    // 19
                                                         % STRING(stats.nLog.F)     // 20
                                                         % STRING(stats.nLog.p)     // 21
                                                         % STRING(stats.nLog.SSM)   // 22
                                                         % STRING(stats.nLog.SSM_D) // 23
                                                         % STRING(stats.nLog.SSE)   // 24
                                                         % STRING(stats.nLog.SSE_D) // 25
                                                         % STRING(stats.nLog.SST)   // 26
                                                         % STRING(stats.nLog.SST_D) // 27
                                                         % STRING(stats.wLog.cor)
                                                         % STRING(stats.wLog.slope)
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

Scripts StatsWriter::inflectSummary(const std::vector<LinearStats> &stats, const Units &units)
{
//    for (const auto &i : stats)
//    {
//        const auto n_lm = i.linear(false);
//        const auto l_lm = i.linear(true);
//        
//        // Calcluate the inflect point after log-transformation
//        const auto inflect = i.inflect(true);
//        
//        // Remember the break-point is on the log-scale, we'll need to convert it back
//        const auto b = pow(2, inflect.b);
//
//        /*
//        return (boost::format(inflectSummary()) % src                          // 1
//                % stats.n_endo
//                % stats.n_chrT
//                % units
//                % stats.hist.size()            // 5
//                % (ref.empty() ? units : ref)  // 6
//                % b
//                % inflect.id
//                % inflect.lInt                 // 9
//                % inflect.lSl                  // 10
//                % inflect.lR2                  // 11
//                % inflect.rInt                 // 12
//                % inflect.rSl                  // 13
//                % inflect.rR2                  // 14
//                % n_lm.r                       // 15
//                % n_lm.m                       // 16
//                % n_lm.r2                      // 17
//                % n_lm.f                       // 18
//                % n_lm.p                       // 19
//                % n_lm.ssm                     // 20
//                % n_lm.ssm_df                  // 21
//                % n_lm.sse                     // 22
//                % n_lm.sse_df                  // 23
//                % n_lm.sst                     // 24
//                % n_lm.sst_df                  // 25
//                % l_lm.r                       // 26
//                % l_lm.m                       // 27
//                % l_lm.r2                      // 28
//                % l_lm.f                       // 29
//                % l_lm.p                       // 30
//                % l_lm.ssm                     // 31
//                % l_lm.ssm_df                  // 32
//                % l_lm.sse                     // 33
//                % l_lm.sse_df                  // 34
//                % l_lm.sst                     // 35
//                % l_lm.sst_df                  // 36
//                ).str();
//        */
//    }
    
    return "";
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