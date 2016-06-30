#include "stats/analyzer.hpp"
#include "writers/r_writer.hpp"

using namespace Anaquin;

// Defined in main.cpp
extern Path __output__;

// Defined in resources.cpp
extern Scripts PlotScatter();

// Defined in resources.cpp
extern Scripts PlotFold();

// Defined in resources.cpp
extern Scripts PlotSensitivity();

// Defined in main.cpp
extern FileName mixture();

// Defined in resources.cpp
extern Scripts PlotVROC();

Scripts RWriter::createVROC(const FileName &file, const std::string &score)
{
    return (boost::format(PlotVROC()) % date()
                                      % __full_command__
                                      % __output__
                                      % file
                                      % score).str();
}

Scripts RWriter::createSensitivity(const FileName    &file,
                                   const std::string &title,
                                   const std::string &xlab,
                                   const std::string &ylab,
                                   const std::string &expected,
                                   const std::string &measured,
                                   bool showLOQ)
{
    return (boost::format(PlotSensitivity()) % date()
                                             % __full_command__
                                             % __output__
                                             % file
                                             % title
                                             % xlab
                                             % ylab
                                             % ("log2(data$" + expected + ")")
                                             % ("data$" + measured)
                                             % (showLOQ ? "TRUE" : "FALSE")).str();
}

Scripts RWriter::createFold(const FileName    &file,
                            const std::string &title,
                            const std::string &xlab,
                            const std::string &ylab,
                            const std::string &expected,
                            const std::string &measured,
                            bool shouldLog)
{
    const auto exp = shouldLog ? ("log2(data$" + expected + ")") : ("data$" + expected);
    const auto obs = shouldLog ? ("log2(data$" + measured + ")") : ("data$" + measured);
    
    return (boost::format(PlotFold()) % date()
                                      % __full_command__
                                      % __output__
                                      % file
                                      % title
                                      % xlab
                                      % ylab
                                      % exp
                                      % obs
                                      % "TRUE").str();
}

Scripts RWriter::createScatterNoLog(const FileName    &file,
                                    const std::string &title,
                                    const std::string &xlab,
                                    const std::string &ylab,
                                    const std::string &expected,
                                    const std::string &measured,
                                    bool showLOQ)
{
    return (boost::format(PlotScatter()) % date()
                                          % __full_command__
                                          % __output__
                                          % file
                                          % title
                                          % xlab
                                          % ylab
                                          % ("data$" + expected)
                                          % ("data$" + measured)
                                          % (showLOQ ? "TRUE" : "FALSE")).str();
}

Scripts RWriter::createMultiScatter(const FileName    &file,
                                    const std::string &title,
                                    const std::string &xlab,
                                    const std::string &ylab,
                                    const std::string &expected,
                                    const std::string &measured,
                                    bool showLOQ,
                                    bool shouldLog)
{
    const auto exp = shouldLog ? ("log2(data$" + expected + ")") : ("data$" + expected);
    const auto obs = shouldLog ? ("log2(data[,2:ncol(data)])") : ("data[,2:ncol(data)]");
    
    return (boost::format(PlotScatter()) % date()
                                         % __full_command__
                                         % __output__
                                         % file
                                         % title
                                         % xlab
                                         % ylab
                                         % exp
                                         % obs
                                         % (showLOQ ? "TRUE" : "FALSE")).str();
}

Scripts RWriter::createScatterNeedLog(const FileName    &file,
                                      const std::string &title,
                                      const std::string &xlab,
                                      const std::string &ylab,
                                      const std::string &expected,
                                      const std::string &measured,
                                      bool showLOQ)
{
    return (boost::format(PlotScatter()) % date()
                                         % __full_command__
                                         % __output__
                                         % file
                                         % title
                                         % xlab
                                         % ylab
                                         % ("log2(data$" + expected + ")")
                                         % ("log2(data$" + measured + ")")
                                         % (showLOQ ? "TRUE" : "FALSE")).str();
}

Scripts RWriter::createScript(const FileName &file, const Scripts &script)
{
    return (boost::format(script) % date()
                                  % __full_command__
                                  % __output__
                                  % file).str();
}

Scripts RWriter::createScript(const FileName &file, const Scripts &script, const std::string &x)
{
    return (boost::format(script) % date()
                                  % __full_command__
                                  % __output__
                                  % file
                                  % x).str();
}

Scripts StatsWriter::linearSummary(const FileName &file,
                                   const FileName &annot,
                                   const LinearStats &stats,
                                   const MappingStats &mStats,
                                   const Hist &hist,
                                   const Units &units)
{
    const auto summary = "   ***\n"
                         "   *** Number of %1% for synthetic and genome detected in the input file\n"
                         "   ***\n\n"
                         "   Synthetic: %2% (%3%)\n"
                         "   Genome:    %4% (%5%)\n\n"
                         "   ***\n"
                         "   *** Reference annotation\n"
                         "   ***\n\n"
                         "   Annotation: %6%\n"
                         "   Mixture:    %7%\n"
                         "   Reference:  %8% %1%\n"
                         "   Detected:   %9% %1%\n\n"
                         "   Limit:      %10% (%11%)\n\n"
                         "   ***\n"
                         "   ***   Correlation: Pearson’s correlation\n"
                         "   ***   Slope:       Regression slope for the regression\n"
                         "   ***   R2:          Coefficient of determination for the regression\n"
                         "   ***   F-stat:      F-test statistic under the null hypothesis\n"
                         "   ***   P-value:     P-value under the null hypothesis\n"
                         "   ***   SSM:         Total sum of squares for model\n"
                         "   ***   SSE:         Total sum of squares for residuals\n"
                         "   ***   SST:         Total sum of squares\n"
                         "   ***\n\n"
                         "   Correlation: %12%\n"
                         "   Slope:       %13%\n"
                         "   R2:          %14%\n"
                         "   F-statistic: %15%\n"
                         "   P-value:     %16%\n"
                         "   SSM:         %17%, DF: %18%\n"
                         "   SSE:         %19%, DF: %20%\n"
                         "   SST:         %21%, DF: %22%\n\n"
                         "   ***\n"
                         "   *** The following statistics are computed on the log2 scale.\n"
                         "   ***\n"
                         "   ***   Eg: If the data points are (1,1), (2,2). The correlation will\n"
                         "   ***       be computed on (log2(1), log2(1)), (log2(2), log2(2)))\n"
                         "   ***\n\n"
                         "   Correlation: %23%\n"
                         "   Slope:       %24%\n"
                         "   R2:          %25%\n"
                         "   F-statistic: %26%\n"
                         "   P-value:     %27%\n"
                         "   SSM:         %28%, DF: %29%\n"
                         "   SSE:         %30%, DF: %31%\n"
                         "   SST:         %32%, DF: %33%\n";
    
    const auto nlm = stats.linear(false);
    const auto llm = stats.linear(true);

    return (boost::format(summary) % units
                                   % mStats.n_syn
                                   % mStats.synProp()
                                   % mStats.n_gen
                                   % mStats.genProp()
                                   % annot
                                   % mixture()           // 7
                                   % hist.size()         // 8
                                   % count(hist)         // 9
                                   % stats.limit.id      // 10
                                   % stats.limit.abund   // 11
                                   % nlm.r               // 12
                                   % nlm.m
                                   % nlm.R2
                                   % nlm.F
                                   % nlm.p
                                   % nlm.SSM
                                   % nlm.SSM_D
                                   % nlm.SSE
                                   % nlm.SSE_D
                                   % nlm.SST
                                   % nlm.SST_D
                                   % llm.r                // 23
                                   % llm.m
                                   % llm.R2
                                   % llm.F
                                   % llm.p
                                   % llm.SSM
                                   % llm.SSM_D
                                   % llm.SSE
                                   % llm.SSE_D
                                   % llm.SST
                                   % llm.SST_D).str();
}

Scripts StatsWriter::inflectSummary()
{
    return "Summary for input: %1%\n\n"
           "%41%"
           "   ***\n"
           "   *** Number of %40% for synthetic and genome detected in the input file\n"
           "   ***\n\n"
           "   Synthetic: %2% (%3%)\n"
           "   Genome:    %4% (%5%)\n\n"
           "   ***\n"
           "   *** Reference annotation\n"
           "   ***\n\n"
           "   Annotation: %6%\n"
           "   Mixture:    %42%\n"
           "   Reference:  %7% %8%\n"
           "   Detected:   %9% %8%\n\n"
           "   ***\n"
           "   ***   Correlation: Pearson’s correlation\n"
           "   ***   Slope:       Regression slope for the regression\n"
           "   ***   R2:          Coefficient of determination for the regression\n"
           "   ***   F-stat:      F-test statistic under the null hypothesis\n"
           "   ***   P-value:     P-value under the null hypothesis\n"
           "   ***   SSM:         Total sum of squares for model\n"
           "   ***   SSE:         Total sum of squares for residuals\n"
           "   ***   SST:         Total sum of squares\n"
           "   ***\n\n"
           "   ***\n"
           "   *** Limit of Quantification (LOQ). Estimated by piecewise segmented regression.\n"
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
           "   *** Linear regression (log2 scale)\n"
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

Scripts StatsWriter::inflectSummary(const FileName &chrTR, const FileName &endoR, const SInflectStats &stats, const Units &units)
{
    const auto intro_1 = "   ***\n"
                         "   *** The statistics are shown by arithmetic average and standard deviation. For example,\n"
                         "   *** 5.12 ± 0.52 has an arithmetic average of 5.12 and standard deviation 0.52.\n"
                         "   ***\n\n";

    std::string intro;

    if (stats.files.size() > 1)
    {
        intro = intro_1;
    }

    throw "Not Implemented";
    
//    return (boost::format(StatsWriter::inflectSummary()) % STRING(stats.files)      // 1
//                                                         % STRING(stats.n_syn)     // 2
//                                                         % STRING(stats.p_chrT)     // 3
//                                                         % STRING(stats.n_gen)     // 4
//                                                         % STRING(stats.p_endo)     // 5
//                                                         % chrTR                    // 6
//                                                         % stats.n_ref              // 7
//                                                         % STRING(stats.units)      // 8
//                                                         % STRING(stats.n_det)      // 9
//                                                         % STRING(stats.b)          // 10
//                                                         % STRING(stats.bID)        // 11
//                                                         % STRING(stats.lInt)       // 12
//                                                         % STRING(stats.lSl)        // 13
//                                                         % STRING(stats.lR2)        // 14
//                                                         % STRING(stats.rInt)       // 15
//                                                         % STRING(stats.rSl)        // 16
//                                                         % STRING(stats.rR2)        // 17
//                                                         % STRING(stats.nLog.r)     // 18
//                                                         % STRING(stats.nLog.sl)    // 19
//                                                         % STRING(stats.nLog.R2)    // 20
//                                                         % STRING(stats.nLog.F)     // 21
//                                                         % STRING(stats.nLog.p)     // 22
//                                                         % STRING(stats.nLog.SSM)   // 23
//                                                         % STRING(stats.nLog.SSM_D) // 24
//                                                         % STRING(stats.nLog.SSE)   // 25
//                                                         % STRING(stats.nLog.SSE_D) // 26
//                                                         % STRING(stats.nLog.SST)   // 27
//                                                         % STRING(stats.nLog.SST_D) // 28
//                                                         % STRING(stats.wLog.r)
//                                                         % STRING(stats.wLog.sl)
//                                                         % STRING(stats.wLog.R2)
//                                                         % STRING(stats.wLog.F)
//                                                         % STRING(stats.wLog.p)
//                                                         % STRING(stats.wLog.SSM)
//                                                         % STRING(stats.wLog.SSM_D)
//                                                         % STRING(stats.wLog.SSE)
//                                                         % STRING(stats.wLog.SSE_D)
//                                                         % STRING(stats.wLog.SST)
//                                                         % STRING(stats.wLog.SST_D)
//                                                         % units                     // 40
//                                                         % intro
//                                                         % mixture()
//            ).str();
};

Scripts StatsWriter::inflectSummary(const FileName    &ref,
                                    const FileName    &query,
                                    const Hist        &hist,
                                    const LinearStats &stats,
                                    const Units &units)
{
    throw "Not Implemented";
}

Scripts StatsWriter::inflectSummary(const FileName     &chrTR,
                                    const FileName     &endoR,
                                    const FileName     &file,
                                    const Hist         &hist,
                                    const MappingStats &mStats,
                                    const LinearStats  &stats,
                                    const Units &units)
{
    return inflectSummary(chrTR,
                          endoR,
                          std::vector<FileName>     { file   },
                          //std::vector<Hist>         { hist   },
                          std::vector<MappingStats> { mStats },
                          std::vector<LinearStats>  { stats  },
                          units);
}

SInflectStats StatsWriter::multiInfect(const FileName                  &chrTR,
                                       const FileName                  &endoR,
                                       const std::vector<FileName>     &files,
                                       //const std::vector<SequinHist>   &hist,
                                       const std::vector<MappingStats> &mStats,
                                       const std::vector<LinearStats>  &lstats)
{
    SInflectStats r;
    
    for (auto i = 0; i < lstats.size(); i++)
    {
        r.files.add(files[i]);
        
        // Linear regression with logarithm
        const auto l_lm = lstats[i].linear(true);
        
        // Calcluate the inflection point with logarithm
        const auto inf = lstats[i].limitQuant(true);
        
        // Remember the break-point is on the log2-scale, we'll need to convert it back
        const auto b = pow(2, inf.b);
        
        r.n_syn.add((unsigned)mStats[i].n_syn);
        r.n_gen.add((unsigned)mStats[i].n_gen);
        r.p_chrT.add(mStats[i].synProp());
        r.p_endo.add(mStats[i].genProp());
        
        //r.n_ref = hist[i].size();
        //r.n_det.add((unsigned)detect(hist[i]));
        
        r.b.add(b);
        r.lr.add(inf.lr);
        r.rr.add(inf.rr);
        r.bID.add(inf.id);
        r.lInt.add(inf.lInt);
        r.rInt.add(inf.rInt);
        r.lSl.add(inf.lSl);
        r.rSl.add(inf.rSl);
        r.lR2.add(inf.lR2);
        r.rR2.add(inf.rR2);

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
    
    return r;
}

Scripts StatsWriter::inflectSummary(const FileName                  &chrTR,
                                    const FileName                  &endoR,
                                    const std::vector<FileName>     &files,
                                    //const std::vector<SequinHist>   &hist,
                                    const std::vector<MappingStats> &mStats,
                                    const std::vector<LinearStats>  &lstats,
                                    const Units &units)
{
    SInflectStats r;
    r.units = units;
    
    for (auto i = 0; i < lstats.size(); i++)
    {
        r.files.add(files[i]);

        // Linear regression with logarithm
        const auto l_lm = lstats[i].linear(true);

        // Calcluate the inflection point with logarithm
        const auto inf = lstats[i].limitQuant(true);

        // Remember the break-point is on the log2-scale, we'll need to convert it back
        const auto b = pow(2, inf.b);

        r.n_syn.add((unsigned)mStats[i].n_syn);
        r.n_gen.add((unsigned)mStats[i].n_gen);
        r.p_chrT.add(mStats[i].synProp());
        r.p_endo.add(mStats[i].genProp());

        //r.n_ref = hist[i].size();
        //r.n_det.add((unsigned)detect(hist[i]));

        r.b.add((unsigned)b);
        r.bID.add(inf.id);
        r.lInt.add(inf.lInt);
        r.rInt.add(inf.rInt);
        r.lSl.add(inf.lSl);
        r.rSl.add(inf.rSl);
        r.lR2.add(inf.lR2);
        r.rR2.add(inf.rR2);

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

    return inflectSummary(chrTR, endoR, r, units);
}