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

// Defined in resources.cpp
extern Scripts PlotScatterPool();

Scripts RWriter::scatterPool(const Path &path, const FileName &file)
{
    std::stringstream ss;
    ss << PlotScatterPool();

    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % path
                                    % file).str();
}

Scripts StatsWriter::inflectSummary()
{
    return "Summary for input: %1%\n\n"
           "   ***\n"
           "   *** The statistics are shown in arithmetic average and standard deviation. For example,\n"
           "   *** 5.12 ± 0.52 has an arithmetic average of 5.12 and standard deviation 0.52.\n"
           "   ***\n\n"
           "   ***\n"
           "   *** Fraction of genes for synthetic and experiment relative to all genes detected in the input file\n"
           "   ***\n\n"
           "   Synthetic:  %2% (%3%)\n"
           "   Experiment: %4% (%5%)\n\n"
           "   ***\n"
           "   *** Reference annotation (Synthetic)\n"
           "   ***\n\n"
           "   File:      %6%\n"
           "   Reference: %7% %8%\n"
           "   Detected:  %9% %8%\n\n"
    
           "   ***\n"
           "   *** Please refer to the online documentation for more details on the regression statistics.\n"
           "   ***\n"
           "   *** Correlation: Pearson’s correlation\n"
           "   *** Slope:       Regression slope for the regression\n"
           "   *** R2:          Coefficient of determination for the regression\n"
           "   *** F-stat:      The F test statistic under the null hypothesis\n"
           "   *** P-value:     The p-value under the null hypothesis\n"
           "   *** SSM:         Sum of squares of model in ANOVA\n"
           "   *** SSE:         Sum of squares of errors in ANOVA\n"
           "   *** SST:         Total sum of squares in ANOVA\n"
           "   ***\n\n"
           "   ***\n"
           "   *** Limit of Quantificiation (LOQ). Estimated by piecewise segmented regression.\n"
           "   ***\n\n"
           ""
           "   Break: %10% (%11%)\n\n"
           "   ***\n"
           "   *** Below LOQ\n"
           "   ***\n\n"
           "   Intercept: %12%\n"
           "   Slope:     %13%\n"
           "   R2:        %14%\n\n"
           "   ***\n"
           "   *** Above LOQ\n"
           "   ***\n\n"
           "   Intercept: %15%\n"
           "   Slope:     %16%\n"
           "   R2:        %17%\n\n"
           "   ***\n"
           "   *** Overall linear regression\n"
           "   ***\n\n"
           "   Correlation: %18%\n"
           "   Slope:       %19%\n"
           "   R2:          %20%\n"
           "   F-statistic: %21%\n"
           "   P-value:     %22%\n"
           "   SSM:         %23%, DF: %24%\n"
           "   SSE:         %25%, DF: %26%\n"
           "   SST:         %27%, DF: %28%\n\n"
           "   ***\n"
           "   *** Overall linear regression (log2 scale)\n"
           "   ***\n\n"
           "   Correlation: %29%\n"
           "   Slope:       %30%\n"
           "   R2:          %31%\n"
           "   F-statistic: %32%\n"
           "   P-value:     %33%\n"
           "   SSM:         %34%, DF: %35%\n"
           "   SSE:         %36%, DF: %37%\n"
           "   SST:         %38%, DF: %39%\n";
}

Scripts StatsWriter::inflectSummary(const FileName &ref, const SInflectStats &stats)
{    
    return (boost::format(StatsWriter::inflectSummary()) % STRING(stats.files)      // 1
                                                         % STRING(stats.n_chrT)     // 2
                                                         % STRING(stats.p_chrT)     // 3
                                                         % STRING(stats.n_endo)     // 4
                                                         % STRING(stats.p_endo)     // 5
                                                         % ref                      // 6
                                                         % STRING(stats.n_ref)      // 7
                                                         % STRING(stats.units)      // 8
                                                         % STRING(stats.n_det)      // 9
                                                         % STRING(stats.b)          // 10
                                                         % STRING(stats.bID)        // 11
                                                         % STRING(stats.lInt)       // 12
                                                         % STRING(stats.lSl)        // 13
                                                         % STRING(stats.lR2)        // 14
                                                         % STRING(stats.rInt)       // 15
                                                         % STRING(stats.rSl)        // 16
                                                         % STRING(stats.rR2)        // 17
                                                         % STRING(stats.nLog.r)     // 18
                                                         % STRING(stats.nLog.sl)    // 19
                                                         % STRING(stats.nLog.R2)    // 20
                                                         % STRING(stats.nLog.F)     // 21
                                                         % STRING(stats.nLog.p)     // 22
                                                         % STRING(stats.nLog.SSM)   // 23
                                                         % STRING(stats.nLog.SSM_D) // 24
                                                         % STRING(stats.nLog.SSE)   // 25
                                                         % STRING(stats.nLog.SSE_D) // 26
                                                         % STRING(stats.nLog.SST)   // 27
                                                         % STRING(stats.nLog.SST_D) // 28
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

Scripts StatsWriter::inflectSummary(const FileName                  &ref,
                                    const std::vector<FileName>     &files,
                                    const std::vector<MappingStats> &mStats,
                                    const std::vector<LinearStats>  &stats,
                                    const Units &units)
{
    SInflectStats r;
    r.units = units;
    
    for (auto i = 0; i < stats.size(); i++)
    {
        r.files.add(files[i]);

        // Linear regression without logarithm
        const auto n_lm = stats[i].linear(false);
        
        // Linear regression with logarithm
        const auto l_lm = stats[i].linear(true);

        // Calcluate the inflect point after log-transformation
        const auto inf = stats[i].limitQuant(true);
        
        // Remember the break-point is on the log2-scale, we'll need to convert it back
        const auto b = pow(2, inf.b);

        r.n_chrT.add(mStats[i].n_chrT);
        r.n_endo.add(mStats[i].n_endo);
        r.p_chrT.add(mStats[i].chrTProp());
        r.p_endo.add(mStats[i].endoProp());

        r.n_ref.add(mStats[i].hist.size());
        r.n_det.add(detect(mStats[i].hist));
        
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
    
    return inflectSummary(ref, r);
}

Scripts RWriter::createSplice(const Path &working, const FileName &fpkms)
{
    std::stringstream ss;
    ss << PlotSplice();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % working
                                    % fpkms).str();
}

Scripts RWriter::createLODR(const Path &working, const FileName &dFile, const std::string &cFile, const std::string &lvl)
{
    std::stringstream ss;
    ss << PlotLODR();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % working
                                    % dFile).str();
}

Scripts RWriter::createMA(const Path &working, const FileName &file, const std::string &lvl)
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