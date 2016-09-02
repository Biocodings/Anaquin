#include <cmath>
#include "tools/gtf_data.hpp"
#include "RnaQuin/r_express.hpp"
#include "parsers/parser_gtf.hpp"
#include "parsers/parser_express.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

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
    if (Standard::isSynthetic(x.cID))
    {
        stats.countSyn++;
        const auto &r = Standard::instance().r_trans;
        
        SequinID id;
        Concent  exp = NAN;
        Measured obs = NAN;

        switch (metrs)
        {
            case Metrics::Isoform:
            {
                const auto m = r.match(x.id);
                
                if (m)
                {
                    if (!isnan(x.abund) && x.abund)
                    {
                        id  = m->id;
                        exp = m->concent(Mix_1);
                        obs = x.abund;
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
                        exp = r.concent(x.id);
                        obs = x.abund;
                    }
                }
                else
                {
                    o.logWarn((boost::format("%1% not found. Unknown sequin gene.") % x.id).str());
                }
                
                break;
            }
        }
        
        if (!id.empty())
        {
            auto &ls = metrs == RExpress::Metrics::Isoform ? stats.isos : stats.genes;
            
            if (ls.contains(id))
            {
                // This happens to Cufflink guided assembly...
                o.warn("Duplicate: " + id);
                
                ls.sum(id, exp, obs);
                stats.countSyn--;
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
        stats.countGen++;

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
            const auto &r = Standard::instance().r_trans;

            if (stats.genes.empty() && !stats.isos.empty())
            {
                std::map<GeneID, FPKM> express;
                
                for (const auto &i : stats.isos)
                {
                    const auto m = r.findTrans(ChrIS, i.first);
                    
                    // Add up the isoform expression for each of the sequin gene
                    express[m->gID] += i.second.y;
                }
                
                // Important, we'll need to reset counting for isoforms
                stats.countSyn = 0;
                
                // Start again at the gene level
                stats.limit = Limit();
                
                for (const auto &i : express)
                {
                    const auto input = r.concent(i.first);
                    
                    stats.countSyn++;
                    stats.genes.add(i.first, input, i.second);
                    
                    if (isnan(stats.limit.abund) || input < stats.limit.abund)
                    {
                        stats.limit.id = i.first;
                        stats.limit.abund = input;
                    }
                }

                o.logInfo("Sequin isoform expressions added for genes: " + std::to_string(stats.genes.size()));
            }
        }
    });
}

static Scripts multipleCSV(const std::vector<RExpress::Stats> &stats, Metrics metrs)
{
    const auto &r = Standard::instance().r_trans;
    
    std::set<SequinID> seqs;
    
    // This is the data structure that will be useful
    std::map<unsigned, std::map<SequinID, Concent>> data;
    
    // Expected concentration
    std::map<SequinID, Concent> expect;
    
    std::stringstream ss;
    ss << "ID\tLength\tInputConcent";
    
    for (auto i = 0; i < stats.size(); i++)
    {
        auto &ls = metrs == RExpress::Metrics::Isoform ? stats[i].isos : stats[i].genes;

        ss << ((boost::format("\tObserved%1%") % (i+1)).str());

        for (const auto &j : ls)
        {
            seqs.insert(j.first);
            expect[j.first]  = j.second.x;
            data[i][j.first] = j.second.y;
            
            A_ASSERT(expect[j.first], "Zero expect concentration");
        }
    }
    
    ss << "\n";
    
    for (const auto &seq : seqs)
    {
        Locus l;
        
        switch (metrs)
        {
            case RExpress::Metrics::Gene:    { l = r.findGene(ChrIS, seq)->l;  break; }
            case RExpress::Metrics::Isoform: { l = r.findTrans(ChrIS, seq)->l; break; }
        }
        
        ss << ((boost::format("%1%\t%2%\t%3%") % seq
                                               % l.length()
                                               % expect.at(seq)).str());
        
        for (auto i = 0; i < stats.size(); i++)
        {
            if (data[i].count(seq))
            {
                ss << "\t" << data[i][seq];
            }
            else
            {
                ss << "\tNA";
            }
        }
        
        ss << "\n";
    }
    
    return ss.str();
}

static void generateSummary(const FileName &summary,
                            const std::vector<FileName >&files,
                            const std::vector<RExpress::Stats> &stats,
                            const RExpress::Options  &o,
                            const Units &units)
{
    const auto &r = Standard::instance().r_trans;
    
    o.info("Generating " + summary);
    o.writer->open(summary);
    
    std::vector<SequinHist>   hists;
    std::vector<LinearStats>  lStats;
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
    
    const auto ms = StatsWriter::multiInfect(files, mStats, lStats);
    
    // Breakpoint estimated by piecewise regression
    const auto b = ms.b.mean();
    
//    // Number of genomic features above the breakpoint
//    SCounts n_above;
//    
//    // Number of genomic features below the breakpoint
//    SCounts n_below;
//    
//    // Counting all replicates...
//    for (const auto &i : stats)
//    {
//        Counts above = 0;
//        Counts below = 0;
//
//        // Counting for each genome gene/isoform...
//        for (const auto &j : i.gData)
//        {
//            assert(!isnan(j.second.abund));
//            
//            if (j.second.abund >= b)
//            {
//                above++;
//            }
//            else
//            {
//                below++;
//            }
//        }
//        
//        n_above.add((Counts)above);
//        n_below.add((Counts)below);
//    }
    
    // No reference coordinate annotation given here
    const auto rSyn = o.metrs == Metrics::Gene || shouldAggregate(o) ? r.countGeneSeqs() : r.countSeqs();

//    const auto hasLOQ = !isnan(ms.b.mean());

    const auto format = "-------RnaExpression Output\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference Transcript Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------%4%\n\n"
                        "       Synthetic: %5%\n"
                        "       Detection Sensitivity: %6% (attomol/ul) (%7%)\n\n"
                        "       Genome: %8%\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %9%\n"
                        "       Correlation: %10%\n"
                        "       R2:          %11%\n"
                        "       F-statistic: %12%\n"
                        "       P-value:     %13%\n"
                        "       SSM:         %14%, DF: %15%\n"
                        "       SSE:         %16%, DF: %17%\n"
                        "       SST:         %18%, DF: %19%\n";
    
    o.writer->write((boost::format(format) % STRING(ms.files)      // 1
                                           % rSyn                  // 2
                                           % MixRef()              // 3
                                           % title                 // 4
                                           % STRING(ms.countSyn)   // 5
                                           % limit.abund           // 6
                                           % limit.id              // 7
                                           % STRING(ms.countGen)   // 8
                                           % STRING(ms.wLog.sl)    // 21
                                           % STRING(ms.wLog.r)     // 22
                                           % STRING(ms.wLog.R2)    // 23
                                           % STRING(ms.wLog.F)     // 24
                                           % STRING(ms.wLog.p)     // 25
                                           % STRING(ms.wLog.SSM)   // 26
                                           % STRING(ms.wLog.SSM_D) // 27
                                           % STRING(ms.wLog.SSE)   // 28
                                           % STRING(ms.wLog.SSE_D) // 29
                                           % STRING(ms.wLog.SST)   // 30
                                           % STRING(ms.wLog.SST_D) // 31
                     ).str());
    o.writer->close();
}

static void generateR(const FileName &output,
                      const FileName &csv,
                      const std::vector<RExpress::Stats> &stats,
                      const RExpress::Options &o)
{
    o.info("Generating " + output);
    o.writer->open(output);
    
    const auto title = o.metrs == Metrics::Gene ? "Gene Expression" : "Isoform Expression";
    
    Units measured;
    
    switch (o.format)
    {
        case RExpress::Format::GTF:      { measured = "FPKM";         break; }
        case RExpress::Format::Text:     { measured = "FPKM";         break; }
        case RExpress::Format::Kallisto: { measured = "K-mer Counts"; break; }
    }

    if (stats.size() == 1)
    {
        o.writer->write(RWriter::createLinear(csv,
                                              title,
                                              "Input Concentration (log2)",
                                              measured + " (log2)",
                                              "InputConcent",
                                              "Observed",
                                              "input",
                                              true));
    }
    else
    {
        o.writer->write(RWriter::createMultiLinear(csv,
                                                   title,
                                                   "Input Concentration (log2)",
                                                   measured + " (log2)",
                                                   "InputConcent",
                                                   "Observed",
                                                   "input",
                                                   true,
                                                   true));
    }
    
    o.writer->close();
}

static void generateCSV(const FileName &output, const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    o.info("Generating " + output);
    o.writer->open(output);

    if (stats.size() == 1)
    {
        const auto format = "%1%\t%2%\t%3%\t%4%";
        
        o.writer->write((boost::format(format) % "ID"
                                               % "Length"
                                               % "InputConcent"
                                               % "Observed").str());
        
        auto &ls = o.metrs == RExpress::Metrics::Isoform ? stats[0].isos : stats[0].genes;
        
        for (const auto &i : ls)
        {
            Locus l;
            
            switch (o.metrs)
            {
                case RExpress::Metrics::Gene:    { l = r.findGene(ChrIS, i.first)->l;  break; }
                case RExpress::Metrics::Isoform: { l = r.findTrans(ChrIS, i.first)->l; break; }
            }

            assert(l.length() > 1);
            o.writer->write((boost::format(format) % i.first
                                                   % l.length()
                                                   % i.second.x
                                                   % i.second.y).str());
        }
    }
    else
    {
        o.writer->write(multipleCSV(stats, o.metrs));
    }

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

    for (const auto &i : stats)
    {
        o.logInfo("Genome: " + toString(i.gData.size()));
    }

    /*
     * Generating RnaExpression_summary.stats
     */
    
    generateSummary("RnaExpression_summary.stats", files, stats, o, units);
    
    /*
     * Generating RnaExpression_sequins.csv
     */
    
    generateCSV("RnaExpression_sequins.csv", stats, o);

    /*
     * Generating RnaExpression_linear.R
     */
    
    generateR("RnaExpression_linear.R", "RnaExpression_sequins.csv", stats, o);
}