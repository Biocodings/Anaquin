#include <cmath>
#include "tools/system.hpp"
#include "tools/gtf_data.hpp"
#include "RnaQuin/RnaQuin.hpp"
#include "RnaQuin/r_express.hpp"
#include "parsers/parser_gtf.hpp"
#include "parsers/parser_express.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

struct MultiStats
{
    SStrings files;
    
    // Eg: 2387648 (15.56%)
    SCounts nSeqs, nEndo;
    
    // Linear regression with logarithm
    SLinearStats stats;
};

// Defined in resources.cpp
extern FileName MixRef();

typedef RExpress::Metrics Metrics;

static bool shouldAggregate(const RExpress::Options &o)
{
    return (o.metrs == Metrics::Gene && o.format == RExpress::Format::GTF) ||
           (o.metrs == Metrics::Gene && o.format == RExpress::Format::Kallisto);
}

template <typename T> void update(RExpress::Stats &stats,
                                  const T &x,
                                  RExpress::Metrics metrs,
                                  const RExpress::Options &o)
{
    const auto &r = Standard::instance().r_rna;
    const auto &l = r.seqsL1();
    
    if (isRNARevChr(x.cID))
    {
        SequinID id;
        Concent  exp = NAN;
        Measured obs = NAN;

        switch (metrs)
        {
            case Metrics::Isoform:
            {
                if (l.count(x.id))
                {
                    if (!isnan(x.abund) && x.abund)
                    {
                        id  = x.id;
                        exp = r.input1(x.id, o.mix);
                        obs = x.abund;
                    }
                    else
                    {
                        o.logWarn((boost::format("Zero or invalid for %1%.") % x.id).str());
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found. Unknown sequin.") % x.id).str());
                }
                
                break;
            }

            case Metrics::Gene:
            {
                const auto m = r.findGene(x.cID, x.id);

                if (m)
                {
                    if (!isnan(x.abund) && x.abund)
                    {
                        id  = x.id;
                        exp = r.concent(x.id, o.mix);
                        obs = x.abund;
                    }
                    else
                    {
                        o.logWarn((boost::format("Zero or invalid for %1%.") % x.id).str());
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found. Unknown sequin gene.") % x.id).str());
                }
                
                break;
            }
        }

        /*
         * A sequin that has zero or invalid measurement is equivalenet being undetected.
         */
        
        if (!isnan(exp))
        {
            stats.nSeqs++;
        }

        if (!id.empty())
        {
            auto &ls = metrs == RExpress::Metrics::Isoform ? stats.isos : stats.genes;
            
            if (ls.count(id))
            {
                // This happens to Cufflink guided assembly...
                o.warn("Duplicate: " + id);
                
                ls.sum(id, exp, obs);
                stats.nSeqs--;
            }
            else
            {
                ls.add(id, exp, obs);
            }

            if (isnan(stats.limit.abund) || exp < stats.limit.abund)
            {
                stats.limit.id = id;
                stats.limit.abund = exp;
            }
        }
    }
    else
    {
        stats.nEndo++;

        // We'll need the information to estimate the numbers below and above the LOQ
        stats.gData[x.id].abund = x.abund;
    }
}

template <typename Functor> RExpress::Stats calculate(const RExpress::Options &o, Functor f)
{
    RExpress::Stats stats;
    
    f(stats);

    return stats;
}

RExpress::Stats RExpress::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    return calculate(o, [&](RExpress::Stats &stats)
    {
        switch (o.format)
        {
            case Format::Kallisto:
            {
                ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &x, const ParserProgress &p)
                {
                    if (p.i && !(p.i % 100000))
                    {
                        o.wait(std::to_string(p.i));
                    }
                    
                    // Should we aggregate because this is at the gene level?
                    if (shouldAggregate(o))
                    {
                        update(stats, x, Metrics::Isoform, o);
                    }
                    else
                    {
                        update(stats, x, o.metrs, o);
                    }
                });
                
                break;
            }
                
            case Format::Text:
            {
                ParserExpress::parse(Reader(file), o.metrs == Metrics::Gene,
                                     [&](const ParserExpress::Data &x, const ParserProgress &p)
                {
                    if (p.i && !(p.i % 100000))
                    {
                        o.wait(std::to_string(p.i));
                    }
                    
                    update(stats, x, o.metrs, o);
                });
                
                break;
            }
                
            case Format::GTF:
            {
                ParserExpress::Data t;

                ParserGTF::parse(file, [&](const ParserGTF::Data &x, const std::string &, const ParserProgress &p)
                {
                    if (p.i && !(p.i % 100000))
                    {
                        o.wait(std::to_string(p.i));
                    }
                    
                    bool matched = false;

                    auto f = [&](Metrics metrs)
                    {
                        if (matched)
                        {
                            t.cID   = x.cID;
                            t.id    = metrs == Metrics::Gene ? x.gID : x.tID;
                            t.abund = x.fpkm;
                            
                            update(stats, t, metrs, o);
                        }
                    };

                    Metrics metrs;
                    
                    switch (metrs = o.metrs)
                    {
                        case Metrics::Isoform:
                        {
                            matched = x.type == RNAFeature::Transcript;
                            break;
                        }
                            
                        case Metrics::Gene:
                        {
                            matched = x.type == RNAFeature::Gene;
                            break;
                        }
                    }
                    
                    f(metrs);

                    /*
                     * Transcriptome GTF file might not have genes defined. We'll need to add
                     * the transcripts together.
                     */

                    if (shouldAggregate(o))
                    {
                        matched = x.type == RNAFeature::Transcript;
                        f(Metrics::Isoform);
                    }
                });

                break;
            }
        }
        
        if (shouldAggregate(o))
        {
            const auto &r = Standard::instance().r_rna;

            if (stats.genes.empty() && !stats.isos.empty())
            {
                std::map<GeneID, FPKM> express;
                
                for (const auto &i : stats.isos)
                {
                    // Add up the isoform expression for the sequin genes
                    express[isoform2Gene(i.first)] += i.second.y;
                }
                
                // Important, we'll need to reset counting for isoforms
                stats.nSeqs = 0;
                
                // Start again at the gene level
                stats.limit = Limit();
                
                for (const auto &i : express)
                {
                    // Input concentration at the gene level
                    const auto input = r.input2(i.first, o.mix);
                    
                    stats.nSeqs++;
                    stats.genes.add(i.first, input, i.second);
                    
                    if (isnan(stats.limit.abund) || input < stats.limit.abund)
                    {
                        stats.limit.id = i.first;
                        stats.limit.abund = input;
                    }
                }
            }
        }
    });
}

static Scripts multipleCSV(const std::vector<RExpress::Stats> &stats, Metrics metrs)
{
    std::set<SequinID> seqs;
    
    // This is the data structure that will be useful
    std::map<unsigned, std::map<SequinID, Concent>> data;
    
    // Expected concentration
    std::map<SequinID, Concent> expect;
    
    std::stringstream ss;
    ss << "ID\tInput";
    
    for (auto i = 0; i < stats.size(); i++)
    {
        auto &ls = metrs == RExpress::Metrics::Isoform ? stats[i].isos : stats[i].genes;

        ss << ((boost::format("\tObserved%1%") % (i+1)).str());

        for (const auto &j : ls)
        {
            seqs.insert(j.first);
            expect[j.first]  = j.second.x;
            data[i][j.first] = j.second.y;
            
            A_CHECK(expect[j.first], "Zero expect concentration");
        }
    }
    
    ss << "\n";
    
    for (const auto &seq : seqs)
    {
        ss << ((boost::format("%1%\t%2%") % seq % expect.at(seq)).str());
        
        for (auto i = 0; i < stats.size(); i++)
        {
            if (data[i].count(seq))
            {
                ss << "\t" << data[i][seq];
            }
            else
            {
                ss << "\tNaN";
            }
        }
        
        ss << "\n";
    }
    
    return ss.str();
}

static MultiStats multiStats(const std::vector<FileName>     &files,
                             const std::vector<MappingStats> &mStats,
                             const std::vector<SequinStats>  &lstats)
{
    A_ASSERT(files.size()  == mStats.size());
    A_ASSERT(mStats.size() == lstats.size());

    MultiStats r;
    
    for (auto i = 0; i < lstats.size(); i++)
    {
        const auto lm = lstats[i].linear();
        
        r.nSeqs.add(mStats[i].nSeqs);
        r.nEndo.add(mStats[i].nEndo);
        
        r.files.add(files[i]);
        r.stats.p.add(lm.p);
        r.stats.r.add(lm.r);
        r.stats.F.add(lm.F);
        r.stats.sl.add(lm.m);
        r.stats.R2.add(lm.R2);
        r.stats.SSM.add(lm.SSM);
        r.stats.SSE.add(lm.SSE);
        r.stats.SST.add(lm.SST);
        r.stats.SSM_D.add(lm.SSM_D);
        r.stats.SSE_D.add(lm.SSE_D);
        r.stats.SST_D.add(lm.SST_D);
    }
    
    return r;
}

Scripts RExpress::generateSummary(const std::vector<FileName> &tmp,
                                  const std::vector<RExpress::Stats> &stats,
                                  const RExpress::Options &o,
                                  const Units &units)
{
    const auto files = path2file(tmp);
    const auto &r = Standard::instance().r_rna;
    
    std::vector<SequinHist>   hists;
    std::vector<SequinStats>  lStats;
    std::vector<MappingStats> mStats;
    
    // Detection limit for the replicates
    Limit limit;
    
    for (auto i = 0; i < files.size(); i++)
    {
        auto &ls = o.metrs == RExpress::Metrics::Isoform ? stats[i].isos : stats[i].genes;
        
        mStats.push_back(stats[i]);
        lStats.push_back(ls);
        
        // Not every replicate is defined...
        if (!stats[i].limit.id.empty())
        {
            if (isnan(limit.abund) || stats[i].limit.abund < limit.abund)
            {
                limit = stats[i].limit;
            }
        }
    }
    
    assert(!isnan(limit.abund) && !limit.id.empty());
    
    const auto title = (o.metrs == Metrics::Gene ? "Genes Expressed" : "Isoform Expressed");
    const auto ms    = multiStats(files, mStats, lStats);
    const auto count = o.metrs == Metrics::Gene || shouldAggregate(o) ? r.seqsL1().size() : r.seqsL2().size();
    
    const auto format = "-------RnaExpression Output\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference Transcript Annotations\n\n"
                        "       Synthetic: %2% %3%\n"
                        "       Mixture file: %4%\n\n"
                        "-------%5%\n\n"
                        "       Sequin: %6%\n"
                        "       Detection Sensitivity: %7% (attomol/ul) (%8%)\n\n"
                        "       Genome: %9%\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %10%\n"
                        "       Correlation: %11%\n"
                        "       R2:          %12%\n"
                        "       F-statistic: %13%\n"
                        "       P-value:     %14%\n"
                        "       SSM:         %15%, DF: %16%\n"
                        "       SSE:         %17%, DF: %18%\n"
                        "       SST:         %19%, DF: %20%\n";
    
    return (boost::format(format) % STRING(ms.files)       // 1
                                  % count                  // 2
                                  % units                  // 3
                                  % MixRef()               // 4
                                  % title                  // 5
                                  % STRING(ms.nSeqs)       // 6
                                  % limit.abund            // 7
                                  % limit.id               // 8
                                  % STRING(ms.nEndo)       // 9
                                  % STRING(ms.stats.sl)    // 10
                                  % STRING(ms.stats.r)     // 11
                                  % STRING(ms.stats.R2)    // 12
                                  % STRING(ms.stats.F)     // 13
                                  % STRING(ms.stats.p)     // 14
                                  % STRING(ms.stats.SSM)   // 15
                                  % STRING(ms.stats.SSM_D) // 16
                                  % STRING(ms.stats.SSE)   // 17
                                  % STRING(ms.stats.SSE_D) // 18
                                  % STRING(ms.stats.SST)   // 19
                                  % STRING(ms.stats.SST_D) // 20
                     ).str();
}

void RExpress::writeSummary(const FileName &file,
                            const std::vector<FileName >&srcs,
                            const std::vector<RExpress::Stats> &stats,
                            const RExpress::Options &o,
                            const Units &units)
{
    o.info("Generating " + file);
    o.writer->open(file);
    o.writer->write(RExpress::generateSummary(srcs, stats, o, units));
    o.writer->close();
}

Scripts RExpress::generateRLinear(const FileName &csv,
                                  const std::vector<RExpress::Stats> &stats,
                                  const RExpress::Options &o)
{
    const auto title = o.metrs == Metrics::Gene ? "Gene Expression" : "Isoform Expression";
    
    Units measured;
    
    switch (o.format)
    {
        case RExpress::Format::GTF:      { measured = "FPKM"; break; }
        case RExpress::Format::Text:     { measured = "FPKM"; break; }
        case RExpress::Format::Kallisto: { measured = "Transcripts per million"; break; }
    }

    if (stats.size() == 1)
    {
        return RWriter::createRLinear(csv,
                                      o.work,
                                      title,
                                      "Input Concentration (log2)",
                                       measured + " (log2)",
                                      "log2(data$InputConcent)",
                                      "log2(data[,3:ncol(data)])",
                                      "input",
                                       true);
    }
    else
    {
        return RWriter::createRLinear(csv,
                                      o.work,
                                      title,
                                      "Input Concentration (log2)",
                                       measured + " (log2)",
                                      "log2(data$InputConcent)",
                                      "log2(data[,3:ncol(data)])",
                                      "input",
                                       true);
    }
}

void RExpress::writeRLinear(const FileName &file,
                            const FileName &csv,
                            const std::vector<RExpress::Stats> &stats,
                            const RExpress::Options &o)
{
    o.info("Generating " + file);
    o.writer->open(file);
    o.writer->write(RExpress::generateRLinear(csv, stats, o));
    o.writer->close();
}

Scripts RExpress::generateCSV(const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    if (stats.size() == 1)
    {
        std::stringstream ss;
        
        const auto format = "%1%\t%2%\t%3%\n";
        ss << (boost::format(format) % "ID" % "Input" % "Observed").str();
        
        auto &x = o.metrs == RExpress::Metrics::Isoform ? stats[0].isos : stats[0].genes;
        
        for (const auto &i : x)
        {
            ss << (boost::format(format) % i.first
                                         % i.second.x
                                         % i.second.y).str();
        }

        return ss.str();
    }
    else
    {
        return multipleCSV(stats, o.metrs);
    }
}

void RExpress::writeCSV(const FileName &output, const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    o.info("Generating " + output);
    o.writer->open(output);
    o.writer->write(generateCSV(stats, o));
    o.writer->close();
}

void RExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto m = std::map<RExpress::Metrics, std::string>
    {
        { RExpress::Metrics::Gene,    "genes"    },
        { RExpress::Metrics::Isoform, "isoforms" },
    };
    
    switch (o.metrs)
    {
        case Metrics::Gene:    { o.info("Gene Expresssion");   break; }
        case Metrics::Isoform: { o.info("Isoform Expression"); break; }
    }
    
    const auto units = m.at(o.metrs);
    const auto stats = analyze(files, o);

    /*
     * Generating RnaExpression_summary.stats
     */
    
    RExpress::writeSummary("RnaExpression_summary.stats", files, stats, o, units);
    
    /*
     * Generating RnaExpression_sequins.csv
     */
    
    RExpress::writeCSV("RnaExpression_sequins.csv", stats, o);

    /*
     * Generating RnaExpression_linear.R
     */
    
    RExpress::writeRLinear("RnaExpression_linear.R", "RnaExpression_sequins.csv", stats, o);
}
