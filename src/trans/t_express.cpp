#include "data/experiment.hpp"
#include "trans/t_express.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_tracking.hpp"
#include "parsers/parser_stringtie.hpp"
#include <ss/regression/segmented.hpp>

using namespace Anaquin;

typedef TExpress::Level    Level;
typedef TExpress::Software Software;

template <typename T> void update(TExpress::Stats &stats, const T &t, const TExpress::Options &o)
{
    if (t.cID != Standard::chrT)
    {
        stats.n_endo++;
    }
    else
    {
        stats.n_chrT++;
        
    }

    if (t.cID == Standard::chrT)
    {
        const auto &r = Standard::instance().r_trans;
        
        switch (o.lvl)
        {
            case Level::Isoform:
            {
                const TransData *m = nullptr;
                
                // Try to match by name if possible
                m = r.match(t.id);
                
                if (!m)
                {
                    // Try to match by locus (de-novo assembly)
                    m = r.match(t.l, Overlap);
                }
                
                if (!m)
                {
                    o.logWarn((boost::format("%1% not found. Unknown isoform.") % t.id).str());
                }
                else
                {
                    stats.hist.at(m->id)++;
                    
                    if (t.fpkm)
                    {
                        stats.data[t.cID].add(t.id, m->abund(Mix_1), t.fpkm);
                    }
                }
                
                break;
            }
                
            case Level::Gene:
            {
                const TransRef::GeneData *m = nullptr;
                
                // Try to match by name if possible
                m = r.findGene(t.cID, t.id);
                
                if (!m)
                {
                    // Try to match by locus (de-novo assembly)
                    m = r.findGene(t.cID, t.l, Contains);
                }
                
                if (m)
                {
                    stats.hist.at(m->id)++;
                    
                    if (t.fpkm)
                    {
                        stats.data[t.cID].add(t.id, m->abund(Mix_1), t.fpkm);
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found. Unknown gene.") % t.id).str());
                }
                
                break;
            }
                
            case Level::Exon:
            {
                throw "Not Implemented";
            }
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

    switch (o.lvl)
    {
        case Level::Exon:
        {
            throw "Not Implemented";
        }

        case Level::Isoform:
        {
            stats.hist = r.hist();
            break;
        }

        case Level::Gene:
        {
            stats.hist = r.geneHist(ChrT);
            break;
        }
    }
    
    f(stats);
    
    switch (o.lvl)
    {
        case Level::Exon:
        {
            throw "Not Implemented";
        }

        case Level::Isoform:
        {
            stats.limit = r.limit(stats.hist);
            break;
        }

        case Level::Gene:
        {
            stats.limit = r.limitGene(stats.hist);
            break;
        }
    }

    return stats;
}

TExpress::Stats TExpress::analyze(const std::vector<Expression> &exps, const Options &o)
{
    return calculate(o, [&](TExpress::Stats &stats)
    {
        for (const auto &i : exps)
        {
            update(stats, i, o);
        }
    });
}

TExpress::Stats TExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);

    return calculate(o, [&](TExpress::Stats &stats)
    {
        switch (o.soft)
        {
            case Software::Cufflinks:
            {
                ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
                {
                    update(stats, t, o);
                });
                
                break;
            }
                
            case Software::StringTie:
            {
                switch (o.lvl)
                {
                    case Level::Isoform:
                    {
                        ParserStringTie::parseIsoforms(file, [&](const ParserStringTie::STExpression &t, const ParserProgress &)
                        {
                            update(stats, t, o);
                        });

                        break;
                    }
                        
                    case Level::Gene:
                    {
                        ParserStringTie::parseGenes(file, [&](const ParserStringTie::STExpression &t, const ParserProgress &)
                        {
                            update(stats, t, o);
                        });

                        break;
                    }
                        
                    case Level::Exon:
                    {
                        throw "Not Implemented";
                    }
                }
                
                break;
            }
        }
    });
}

static void writeSummary(const TExpress::Stats   &stats,
                         const FileName          &file,
                         const std::string       &name,
                         const ChromoID          &cID,
                         const std::string       &units,
                         const TExpress::Options &o)
{
    const auto sample = extractFile(file);
    
    // Create the directory if haven't
    o.writer->create(name);
    
    o.writer->open(name + "/TransExpress_summary.stats");
    o.writer->write(StatsWriter::linearInflect(file, stats, cID, units));
    o.writer->close();
}

static void writeScatter(const TExpress::Stats   &stats,
                         const FileName          &file,
                         const std::string       &name,
                         const ChromoID          &cID,
                         const std::string       &units,
                         const TExpress::Options &o)
{
    const auto sample = extractFile(file);
    
    // Create the directory if haven't
    o.writer->create(name);

    o.writer->open(name + "/TransExpress_scatter.R");
    
    o.writer->write(RWriter::scatter(stats,
                                     ChrT,
                                     "",
                                     "TransExpress",
                                     "Expected concentration (attomol/ul)",
                                     "Measured coverage (FPKM)",
                                     "Expected concentration (log2 attomol/ul)",
                                     "Measured coverage (log2 FPKM)"));
    o.writer->close();
}

void TExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TExpress::analyze(files, o);
    
    const auto m = std::map<TExpress::Level, std::string>
    {
        { TExpress::Level::Gene, "gene"    },
        { TExpress::Level::Gene, "exon"    },
        { TExpress::Level::Gene, "isoform" },
    };
    
    const auto units = m.at(o.lvl);

    o.info("Generating statistics");
    
    std::map<ChromoID, Accumulator<double>> accs;    
    
    /*
     * We'll pool the information with accumulators to generate a summary for all replicates.
     */
    
    for (auto i = 0; i < files.size(); i++)
    {
        /*
         * Generating summary statistics for the replicate
         */
        
        writeSummary(stats[i], files[i], o.exp->names().at(i), ChrT, units, o);
        
        /*
         * Generating scatter plot for the replicate
         */

        writeScatter(stats[i], files[i], o.exp->names().at(i), ChrT, units, o);
    }
}