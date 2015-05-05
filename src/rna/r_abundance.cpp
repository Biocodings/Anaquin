#include <iostream>
#include <ss/c.hpp>
#include "expression.hpp"
#include "r_abundance.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_tmap.hpp"
#include "parsers/parser_tracking.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace SS::R;
using namespace Spike;

static bool suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

RAbundanceStats RAbundance::analyze(const std::string &file, const Options &options)
{
    RAbundanceStats stats;
    const auto &s = Standard::instance();

    auto c = RAnalyzer::isoformCounter();

    std::vector<double> x, y;
    std::vector<SequinID> z;

    if (suffix(file, ".tmap"))
    {
        ParserTMap::parse(file, [&](const TMap &t, const ParserProgress &)
        {
            if (s.r_seqs_iA.count(t.refID))
            {
                c[t.refID]++;
                
                if (t.fpkm)
                {
                    const auto &i = s.r_seqs_iA.at(t.refID);

                    x.push_back(i.abund(true));
                    y.push_back(t.fpkm);
                    z.push_back(t.refID);
                }
            }
        });

        stats.s = Expression::analyze(c, s.r_sequin(options.mix));
    }
    else
    {
//        ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &)
//                              {
//                                  assert(s.r_seqs_gA.count(t.geneID));
//                                  
//                                  switch (options.level)
//                                  {
//                                      case Gene:
//                                      {
//                                          c[t.geneID]++;
//                                          const auto &m = s.r_seqs_gA.at(t.geneID);
//                                          
//                                          if (t.fpkm)
//                                          {
//                                              /*
//                                               * The x-axis would be the known concentration for each gene,
//                                               * the y-axis would be the expression (RPKM) reported.
//                                               */
//                                              
//                                              x.push_back(m.abund(true));
//                                              y.push_back(t.fpkm);
//                                          }
//                                          
//                                          break;
//                                      }
//                                          
//                                      case Isoform:
//                                      {
//                                          c[t.trackID]++;
//                                          assert(s.r_seqs_iA.count(t.trackID));
//                                          
//                                          if (t.fpkm)
//                                          {
//                                              const auto &i = s.r_seqs_iA.at(t.trackID);
//                                              
//                                              x.push_back(i.abund(true));
//                                              y.push_back(t.fpkm);                    
//                                          }
//                                          
//                                          break;
//                                      }
//                                  }
//                              });
//        
//        if (options.level == Gene)
//        {
//            stats.s = Expression::analyze(c, s.r_pair(options.mix));
//        }
//        else
//        {
//            stats.s = Expression::analyze(c, s.r_sequin(options.mix));
//        }
    }
    
    assert(!x.empty() && x.size() == y.size() && y.size() == z.size());

    // Perform a linear-model to the abundance
    const auto m = lm("y ~ x", data.frame(SS::c(y), SS::c(x)));

    /*
     * In our analysis, the dependent variable is expression while the independent
     * variable is the known concentraion.
     *
     *     expression = constant + slope * concentraion
     */

    stats.lm.r2 = m.ar2;
    stats.lm.r  = cor(x, y);
    stats.lm.c  = m.coeffs[0].v;
    stats.lm.m  = m.coeffs[1].v;

    const std::string format = "%1%\t%2%\t%3%";
    
    options.writer->open("abundance.stats");
    options.writer->write((boost::format(format) % "r" % "s" % "ss").str());
    options.writer->write((boost::format(format) % stats.lm.r2
                                                 % stats.lm.m
                                                 % stats.s.abund).str());
    options.writer->close();

    /*
     * Generate a plot for the fold-change relationship
     */

    options.writer->open("abundance.R");
    options.writer->write(RWriter::write(x, y, z));
    options.writer->close();

    return stats;
}