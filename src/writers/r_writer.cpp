#include <iostream>
#include "data/experiment.hpp"
#include "writers/r_writer.hpp"

using namespace Anaquin;

Scripts RWriter::createLODR(const std::string &working, const FileName &dFile, const std::string &cFile, const std::string &lvl)
{
    std::stringstream ss;
    ss << PlotLODR();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % working
                                    % cFile
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

Scripts RWriter::createROC(const std::vector<std::string> &seqs, const std::vector<double> &qs)
{
    assert(!seqs.empty() && seqs.size() == qs.size());
    
    std::stringstream ss;
    ss << PlotROC();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % ("\'" + boost::algorithm::join(seqs, "\',\'") + "\'")
                                    % RWriter::concat(qs)
            ).str();
}
