#include "trans/t_express.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_tracking.hpp"
#include "parsers/parser_stringtie.hpp"

using namespace SS;
using namespace Anaquin;

typedef TExpress::Metrics  Metrics;
typedef TExpress::Software Software;

template <typename T> void update(TExpress::Stats &stats, const T &t, const GenericID &id, const TExpress::Options &o)
{
    if (t.cID != Standard::chrT)
    {
        stats.n_expT++;
    }
    else
    {
        stats.n_chrT++;
    }

    const auto &r = Standard::instance().r_trans;

    switch (o.metrs)
    {
        case Metrics::Exon:
        {
            break;
        }
            
        case Metrics::Isoform:
        {
            const TransData *m = nullptr;
            
            // Try to match by name if possible
            m = r.match(id);
            
            if (!m)
            {
                // Try to match by locus (de-novo assembly)
                m = r.match(t.l, Overlap);
            }
            
            if (!m)
            {
                o.logWarn((boost::format("%1% not found. Unknown isoform.") % id).str());
            }
            else
            {
                stats.hist.at(m->id)++;
                
                if (t.fpkm)
                {
                    stats.data[t.cID].add(id, m->abund(Mix_1), t.fpkm);
                }
            }

            break;
        }

        case Metrics::Gene:
        {
            const TransRef::GeneData *m = nullptr;
            
            // Try to match by name if possible
            m = r.findGene(t.cID, id);
            
            if (!m)
            {
                // Try to match by locus (de-novo assembly)
                m = r.findGene(t.cID, t.l, Contains);
            }
            
            if (!m)
            {
                o.logWarn((boost::format("%1% not found. Unknown gene.") % id).str());
            }
            else
            {
                stats.hist.at(m->id)++;
                
                if (t.fpkm)
                {
                    stats.data[t.cID].add(id, m->abund(Mix_1), t.fpkm);
                }
            }

            break;
        }
    }
}

template <typename Functor> TExpress::Stats calculate(const TExpress::Options &o, Functor f)
{
    TExpress::Stats stats;
    
    const auto &r   = Standard::instance().r_trans;
    const auto cIDs = r.chromoIDs();
    
    std::for_each(cIDs.begin(), cIDs.end(), [&](const ChromoID &cID)
    {
        stats.data[cID];
    });

    stats.hist  =  o.metrs == Metrics::Isoform  ? r.hist() : r.geneHist(ChrT);
    stats.limit = (o.metrs == Metrics::Isoform) ? r.limit(stats.hist) : r.limitGene(stats.hist);
    
    f(stats);
    
    return stats;
}

TExpress::Stats TExpress::analyze(const std::vector<Expression> &exps, const Options &o)
{
    return calculate(o, [&](TExpress::Stats &stats)
    {
        for (const auto &i : exps)
        {
            update(stats, i, i.id, o);
        }
    });
}

TExpress::Stats TExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);

    return calculate(o, [&](TExpress::Stats &stats)
    {
        const auto isIsoform = o.metrs == Metrics::Isoform;
        
        switch (o.soft)
        {
            case Software::Cufflinks:
            {
                ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
                {
                    update(stats, t, isIsoform ? t.trackID : t.id, o);
                });
                
                break;
            }
                
            case Software::StringTie:
            {
                switch (o.metrs)
                {
                    case Metrics::Isoform:
                    {
                        ParserStringTie::parseIsoforms(file, [&](const ParserStringTie::STExpression &t, const ParserProgress &)
                        {
                            update(stats, t, t.id, o);
                        });

                        break;
                    }
                        
                    case Metrics::Gene:
                    {
                        ParserStringTie::parseGenes(file, [&](const ParserStringTie::STExpression &t, const ParserProgress &)
                        {
                            update(stats, t, t.id, o);
                        });

                        break;
                    }
                }
                
                break;
            }
        }
    });
}

void TExpress::report(const FileName &file, const Options &o)
{
    const auto stats = TExpress::analyze(file, o);
    const auto units = (o.metrs == Metrics::Isoform) ? "isoforms" : "genes";
    
    o.info("Generating statistics");
    
    for (const auto &i : stats.data)
    {
        /*
         * Generating summary statistics
         */

        o.writer->open("TransExpress_summary.stats");
        o.writer->write(StatsWriter::linear(file, stats, i.first, units));
        o.writer->close();
        
        /*
         * Generating scatter plot
         */
        
        o.writer->open("TransExpress_scatter.R");
        o.writer->write(RWriter::scatter(stats, i.first, "", "TransExpress", "Expected concentration (attomol/ul)", "Measured coverage (FPKM)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 FPKM)"));
        o.writer->close();
    }
}

void TExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TExpress::analyze(files, o);
    const auto units = (o.metrs == Metrics::Isoform) ? "isoforms" : "genes";
    
    /*
     * Generating summary statistics for each replicate
     */
    
    o.info("Generating summary statistics");
    
    for (auto i = 0; i < files.size(); i++)
    {
        //AnalyzeReporter::linear((boost::format("TransExpress_%1%_summary.stats") % files[i]).str(),
          //                      files[i],
            //                    stats[i],
              //                  units,
                //                o.writer);
    }

    /*
     * Generating scatter plot for each replicate
     */
    
    o.info("Generating scatter plot");

    for (auto i = 0; i < files.size(); i++)
    {
        const auto file = (boost::format("TransExpress_%1%") % files[i]).str();
        //AnalyzeReporter::scatter(stats[i], "", file, "Expected concentration (attomol/ul)", "Measured coverage (FPKM)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 FPKM)", o.writer);
    }
    
    /*
     * Generating summary statistics
     */
    
}