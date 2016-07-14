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
                                    const std::string &xname,
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
                                          % xname
                                          % (showLOQ ? "TRUE" : "FALSE")).str();
}

Scripts RWriter::createMultiScatter(const FileName    &file,
                                    const std::string &title,
                                    const std::string &xlab,
                                    const std::string &ylab,
                                    const std::string &expected,
                                    const std::string &measured,
                                    const std::string &xname,
                                    bool showLOQ,
                                    bool shouldLog)
{
    const auto exp = shouldLog ? ("log2(data$" + expected + ")") : ("data$" + expected);
    const auto obs = shouldLog ? ("log2(data[,3:ncol(data)])") : ("data[,3:ncol(data)]");
    
    return (boost::format(PlotScatter()) % date()
                                         % __full_command__
                                         % __output__
                                         % file
                                         % title
                                         % xlab
                                         % ylab
                                         % exp
                                         % obs
                                         % xname
                                         % (showLOQ ? "TRUE" : "FALSE")).str();
}

Scripts RWriter::createScatterNeedLog(const FileName    &file,
                                      const std::string &title,
                                      const std::string &xlab,
                                      const std::string &ylab,
                                      const std::string &expected,
                                      const std::string &measured,
                                      const std::string &xname,
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
                                         % xname
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

SInflectStats StatsWriter::multiInfect(const FileName                  &chrTR,
                                       const FileName                  &endoR,
                                       const std::vector<FileName>     &files,
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