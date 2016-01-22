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

Scripts RWriter::createROC(const std::vector<std::string> &seqs, const std::vector<double> &ps, const std::string &lvl)
{
    assert(!seqs.empty() && seqs.size() == ps.size());

    std::stringstream ss;
    ss << PlotROC();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % ("\'" + boost::algorithm::join(seqs, "\',\'") + "\'")
                                    % RWriter::concat(ps)
                                    % lvl).str();
}
