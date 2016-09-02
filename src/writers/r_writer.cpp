#include "stats/analyzer.hpp"
#include "writers/r_writer.hpp"

using namespace Anaquin;

// Defined in main.cpp
extern Path __output__;

// Defined in resources.cpp
extern Scripts PlotLinear();

// Defined in resources.cpp
extern Scripts PlotFold();

// Defined in resources.cpp
extern Scripts PlotLogistic();

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
    return (boost::format(PlotLogistic()) % date()
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
                            bool shouldLog,
                            const std::string &extra)
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
                                      % "TRUE"
                                      % extra).str();
}

Scripts RWriter::createMultiLinear(const FileName    &file,
                                   const std::string &title,
                                   const std::string &xlab,
                                   const std::string &ylab,
                                   const std::string &expected,
                                   const std::string &measured,
                                   const std::string &xname,
                                   bool showLOQ,
                                   bool shouldLog,
                                   const std::string &extra)
{
    const auto exp = shouldLog ? ("log2(data$" + expected + ")") : ("data$" + expected);
    const auto obs = shouldLog ? ("log2(data[,3:ncol(data)])") : ("data[,3:ncol(data)]");
    
    return (boost::format(PlotLinear()) % date()
                                         % __full_command__
                                         % __output__
                                         % file
                                         % title
                                         % xlab
                                         % ylab
                                         % exp
                                         % obs
                                         % xname
                                         % (showLOQ ? "TRUE" : "FALSE")
                                         % extra).str();
}

Scripts RWriter::createLinear(const FileName    &file,
                              const std::string &title,
                              const std::string &xlab,
                              const std::string &ylab,
                              const std::string &expected,
                              const std::string &measured,
                              const std::string &xname,
                              bool showLOQ,
                              const std::string &extra)
{
    return (boost::format(PlotLinear()) % date()
                                         % __full_command__
                                         % __output__
                                         % file
                                         % title
                                         % xlab
                                         % ylab
                                         % ("log2(data$" + expected + ")")
                                         % ("log2(data$" + measured + ")")
                                         % xname
                                         % (showLOQ ? "TRUE" : "FALSE")
                                         % extra).str();
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

SInflectStats StatsWriter::multiInfect(const std::vector<FileName>     &files,
                                       const std::vector<MappingStats> &mStats,
                                       const std::vector<LinearStats>  &lstats)
{
    SInflectStats r;
    
    for (auto i = 0; i < lstats.size(); i++)
    {
        r.files.add(files[i]);
        
        // Linear regression with logarithm
        const auto l_lm = lstats[i].linear(true);
        
        // Remember the break-point is on the log2-scale, we'll need to convert it back
//        const auto b = pow(2, inf.b);
        
        r.countSyn.add(mStats[i].countSyn);
        r.countGen.add(mStats[i].countGen);
        
//        r.b.add(b);
//        r.lr.add(inf.lr);
//        r.rr.add(inf.rr);
//        r.bID.add(inf.id);
//        r.lInt.add(inf.lInt);
//        r.rInt.add(inf.rInt);
//        r.lSl.add(inf.lSl);
//        r.rSl.add(inf.rSl);
//        r.lR2.add(inf.lR2);
//        r.rR2.add(inf.rR2);

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