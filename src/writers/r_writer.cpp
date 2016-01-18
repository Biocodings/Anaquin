#include <iostream>
#include "data/experiment.hpp"
#include "writers/r_writer.hpp"

using namespace Anaquin;

Scripts RWriter::createLODR(const std::vector<std::string> &seqs,
                            const std::vector<double> &avgs,
                            const std::vector<double> &pvals,
                            const std::vector<double> &logFCs)
{
    std::stringstream ss;
    ss << PlotLODR();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % ("\'" + boost::algorithm::join(seqs, "\',\'") + "\'")
                                    % RWriter::concat(avgs)
                                    % RWriter::concat(pvals)
                                    % RWriter::concat(logFCs)
            ).str();
    
    return ss.str();
}

Scripts RWriter::createMA(const FileName &file, const std::string &lvl)
{
    std::stringstream ss;
    ss << PlotMA();
    
    return (boost::format(ss.str()) % date()
                                    % __full_command__
                                    % file
                                    % lvl).str();
}

//Scripts RWriter::createMA(const CountTable &c)
//{
//    std::stringstream s1;
//    
//    const auto &names = c.names();
//    
//    /*
//     * Construct for the replicates
//     */
//    
//    for (auto i = 0; i < names.size(); i++)
//    {
//        const auto &counts = c.counts(names[i]);
//
//        // Eg: A1 = c(...)
//        s1 << names[i] << "= c(" << RWriter::concat(counts, x2str) << ")";
//
//        if (i != names.size() - 1)
//        {
//            s1 << ",\n";
//        }
//    }
//
//    std::stringstream s2;
//    s2 << PlotMA();
//    
//    return (boost::format(s2.str()) % date()
//                                    % __full_command__
//                                    % s1.str()
//                                    % RWriter::concat(c.ids())
//            ).str();
//}

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
