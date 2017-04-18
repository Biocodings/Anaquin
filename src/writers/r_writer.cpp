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

Scripts RWriter::createLogistic(const FileName    &file,
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
                            const Path        &path,
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
                                      % path
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
                                   const Path        &path,
                                   const std::string &title,
                                   const std::string &xlab,
                                   const std::string &ylab,
                                   const std::string &expected,
                                   const std::string &measured,
                                   const std::string &xname,
                                   bool  showLOQ,
                                   bool  shouldLog,
                                   const std::string &extra)
{
    const auto exp = shouldLog ? ("log2(data$" + expected + ")") : ("data$" + expected);
    const auto obs = shouldLog ? ("log2(data[,3:ncol(data)])") : ("data[,3:ncol(data)]");
    
    return (boost::format(PlotLinear()) % date()
                                         % __full_command__
                                         % path
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

Scripts RWriter::createRLinear(const FileName    &file,
                               const Path        &path,
                               const std::string &title,
                               const std::string &xlab,
                               const std::string &ylab,
                               const std::string &expected,
                               const std::string &measured,
                               const std::string &xname,
                               bool  showLOQ,
                               const std::string &script)
{
    return (boost::format(script.empty() ? PlotLinear() : script)
                                  % date()
                                  % __full_command__
                                  % path
                                  % file
                                  % title
                                  % xlab
                                  % ylab
                                  % expected
                                  % measured
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
