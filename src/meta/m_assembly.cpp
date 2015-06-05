#include "meta/m_assembly.hpp"
#include "parsers/parser_blat.hpp"

using namespace Spike;

MAssembly::Stats MAssembly::analyze(const std::string &file, const Options &options)
{
    MAssembly::Stats stats = Velvet::parse<MAssembly::Stats, MAssembly::AlignedContig>(file);
    
    if (!options.blast.empty())
    {
        std::map<std::string, ParserBlat::BlatLine> psl;

        // Perform a pairwise alignment with blat
        ParserBlat::parse(options.blast, [&](const ParserBlat::BlatLine &l, const ParserProgress &p)
        {
            psl[l.qName] = l;
        });

        std::cout << "Using an aligment file: "  << options.blast << std::endl;
        std::cout << "Creating a abundance plot" << std::endl;

        /*
         * Plot the coverage realtive to the known concentration (in attamoles/ul) of each assembled contig.
         */

        std::vector<double> x, y;
        std::vector<std::string> z;

        for (auto &contig : stats.contigs)
        {
            if (psl.count(contig.id))
            {
                contig.sequin = psl[contig.id].tName;
            }

            if (!contig.sequin.empty())
            {
                //const auto &seq = Standard::instance().m_seq_A.at(node.sequin);
                
                // Known concentration from the matching sequin
                //const auto known = seq.abund();
                
                // Measured coverage
                //const auto measured = node.cov;
                
                //x.push_back(log(known));
                //y.push_back(log(measured));
                //z.push_back(node.sequin);
            }
        }
        
        // Generate a R script for a plot of abundance
        AnalyzeReporter::script("abundance.R", x, y, z, options.writer);
    }
    
    /*
     * Write out assembly results
     */

    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";

    options.writer->open("assembly.stats");
    options.writer->write((boost::format(format) % "Nodes" % "N20" % "N50" % "N80" % "min" % "mean" % "max" % "total" % "reads").str());
    options.writer->write((boost::format(format) % stats.contigs.size()
                                                 % stats.N20
                                                 % stats.N50
                                                 % stats.N80
                                                 % stats.min
                                                 % stats.mean
                                                 % stats.max
                                                 % stats.total
                                                 % -1).str());
    options.writer->close();

    return stats;;
}