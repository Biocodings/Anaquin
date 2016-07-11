#include <stdexcept>
#include "tools/gtf_data.hpp"
#include "RnaQuin/r_express.hpp"
#include "parsers/parser_gtf.hpp"
#include "parsers/parser_express.hpp"

using namespace Anaquin;

typedef RExpress::Metrics Metrics;

template <typename T> void update(RExpress::Stats &stats,
                                  const T &x,
                                  RExpress::Metrics metrs,
                                  const RExpress::Options &o)
{
    if (Standard::isSynthetic(x.cID))
    {
        stats.n_syn++;
        const auto &r = Standard::instance().r_trans;
        
        SequinID id;
        Concent  exp = NAN;
        Measured obs = NAN;

        auto hist = metrs == Metrics::Isoform ? stats.isosHist : stats.geneHist;
        
        switch (metrs)
        {
            case Metrics::Isoform:
            {
                const auto m = r.match(x.id);
                
                if (m)
                {
                    hist.at(x.cID).at(m->id)++;

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
                    hist.at(x.cID).at(x.id)++;
                    
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
            
            ls.add(id, exp, obs);

            if (isnan(ls.limit.abund) || exp < ls.limit.abund)
            {
                ls.limit.id = id;
                ls.limit.abund = exp;
            }
        }
    }
    else
    {
        stats.n_gen++;

        // We'll need the information to estimate the numbers below and above the LOQ
        stats.gData[x.id].abund = x.abund;
    }
}

template <typename Functor> RExpress::Stats calculate(const RExpress::Options &o, Functor f)
{
    RExpress::Stats stats;
    
    const auto &r = Standard::instance().r_trans;
    
    stats.isosHist = r.histIsof();
    stats.geneHist = r.histGene();
    
    assert(!stats.isosHist.empty());
    assert(!stats.geneHist.empty());
    
    f(stats);
    
    if (stats.genes.empty() && stats.isos.empty())
    {
        throw std::runtime_error("Failed to find anything for the synthetic chromosome");
    }
    
    return stats;
}

RExpress::Stats RExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);
    
    return calculate(o, [&](RExpress::Stats &stats)
    {
        switch (o.inputs)
        {
            case Inputs::Text:
            {
                ParserExpress::parse(Reader(file), [&](const ParserExpress::Data &x, const ParserProgress &)
                {
                    update(stats, x, o.metrs, o);
                });
                
                break;
            }
                
            case Inputs::GTF:
            {
                ParserExpress::Data t;

                ParserGTF::parse(file, [&](const ParserGTF::Data &x, const std::string &, const ParserProgress &)
                {
                    bool matched = false;

                    auto f = [&](Metrics metrs)
                    {
                        if (matched)
                        {
                            t.l     = x.l;
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
                     * A transcriptome GTF file might not have genes defined. We'll need to add
                     * the transcripts together.
                     */

                    if (metrs == Metrics::Gene && o.inputs == RExpress::Inputs::GTF)
                    {
                        matched = x.type == RNAFeature::Transcript;
                        f(Metrics::Isoform);
                    }
                });

                break;
            }
        }
        
        if (o.metrs == Metrics::Gene && o.inputs == RExpress::Inputs::GTF)
        {
            const auto &r = Standard::instance().r_trans;

            if (stats.genes.empty() && !stats.isos.empty())
            {
                std::map<GeneID, FPKM> fpkm;
                
                for (const auto &i : stats.isos)
                {
                    const auto m = r.findTrans(ChrIS, i.first);
                    
                    // Add up the isoform expressions for each of the sequin gene
                    fpkm[m->gID] += i.second.y;
                }
                
                for (const auto &i : fpkm)
                {
                    stats.genes.add(i.first, r.concent(i.first), i.second);
                }

                o.info("Sequin isoform expressions added for genes: " + std::to_string(stats.genes.size()));
            }
        }
    });
}

static void writeQueries(const FileName &output, const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{}

static Scripts multipleCSV(const std::vector<RExpress::Stats> &stats, Metrics metrs)
{
    const auto &r = Standard::instance().r_trans;
    
    std::set<SequinID> seqs;
    
    // This is the data structure that will be convenient
    std::map<unsigned, std::map<SequinID, Concent>> data;
    
    // Expected concentration
    std::map<SequinID, Concent> expect;
    
    std::stringstream ss;
    ss << "ID\tLength\tExpected";
    
    for (auto i = 0; i < stats.size(); i++)
    {
        auto &ls = metrs == RExpress::Metrics::Isoform ? stats[i].isos : stats[i].genes;

        ss << ((boost::format("\tObserved%1%") % (i+1)).str());
        
        for (const auto &j : ls)
        {
            seqs.insert(j.first);
            expect[j.first]  = j.second.x;
            data[i][j.first] = j.second.y;
        }
    }
    
    ss << "\n";
    
    for (const auto &seq : seqs)
    {
        ss << ((boost::format("%1%\t%2%\t%3%") % seq
                % r.match(seq)->l.length()
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

/*
 * Generate summary statistics for a single sample and multiple samples.
 */

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
        //hists.push_back(stats[i].hist);
        
        if (isnan(limit.abund) || ls.limit.abund < limit.abund)
        {
            limit = ls.limit;
        }
    }
    
    const auto title = (o.metrs == Metrics::Gene ? "Genes Expressed" : "Isoform Expressed");
    
    const auto ms = StatsWriter::multiInfect(o.rAnnot, o.rAnnot, files, mStats, lStats);
    
    // Breakpoint estimated by piecewise regression
    const auto b = ms.b.mean();
    
    // Number of genomic features above the breakpoint
    SCounts n_above;
    
    // Number of genomic features below the breakpoint
    SCounts n_below;
    
    // Counting all replicates
    for (const auto &i : stats)
    {
        Counts above = 0;
        Counts below = 0;
        
        for (const auto &j : i.gData)
        {
            assert(!isnan(j.second.abund));
            
            if (j.second.abund >= b)
            {
                above++;
            }
            else
            {
                below++;
            }
        }
        
        n_above.add((Counts)above);
        n_below.add((Counts)below);
    }
    
    // No reference coordinate annotation given here
    const auto n_syn = o.metrs == Metrics::Gene ? r.countGeneSeqs() : r.countSeqs();
    
    const auto format = "-------RnaExpression Output\n"
                        "       Summary for input: %1%\n"
                        "       *Arithmetic average and standard deviation are shown\n\n"
                        "-------Reference Transcript Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------%4%\n\n"
                        "       Synthetic: %5%\n"
                        "       Detection Sensitivity: %6% (attomol/ul) (%7%)\n\n"
                        "       Genome: %8%\n\n"
                        "-------Limit of Quantification (LOQ)\n\n"
                        "       *Estimated by piecewise segmented regression\n\n"
                        "       Break LOQ: %9% attomol/ul (%10%)\n\n"
                        "       *Below LOQ\n"
                        "       Intercept:   %11%\n"
                        "       Slope:       %12%\n"
                        "       Correlation: %13%\n"
                        "       R2:          %14%\n"
                        "       Genome:      %15%\n\n"
                        "       *Above LOQ\n"
                        "       Intercept:   %16%\n"
                        "       Slope:       %17%\n"
                        "       Correlation: %18%\n"
                        "       R2:          %19%\n"
                        "       Genome:      %20%\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %21%\n"
                        "       Correlation: %22%\n"
                        "       R2:          %23%\n"
                        "       F-statistic: %24%\n"
                        "       P-value:     %25%\n"
                        "       SSM:         %26%, DF: %27%\n"
                        "       SSE:         %28%, DF: %29%\n"
                        "       SST:         %30%, DF: %31%\n";
    
    o.writer->write((boost::format(format) % STRING(ms.files)      // 1
                                           % n_syn                 // 2
                                           % MixRef()              // 3
                                           % title                 // 4
                                           % STRING(ms.n_syn)      // 5
                                           % limit.abund           // 6
                                           % limit.id              // 7
                                           % STRING(ms.n_gen)      // 8
                                           % STRING(ms.b)          // 9
                                           % STRING(ms.bID)        // 10
                                           % STRING(ms.lInt)       // 11
                                           % STRING(ms.lSl)        // 12
                                           % STRING(ms.lr)         // 13
                                           % STRING(ms.lR2)        // 14
                                           % STRING(n_below)       // 15
                                           % STRING(ms.rInt)       // 16
                                           % STRING(ms.rSl)        // 17
                                           % STRING(ms.rr)         // 18
                                           % STRING(ms.rR2)        // 19
                                           % STRING(n_above)       // 20
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

static void generateCSV(const FileName &output, const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    o.info("Generating " + output);
    o.writer->open(output);
    
    auto &ls = o.metrs == RExpress::Metrics::Isoform ? stats[0].isos : stats[0].genes;
    
    if (stats.size() == 1)
    {
        const auto format = "%1%\t%2%\t%3%\t%4%";
        
        o.writer->write((boost::format(format) % "ID"
                                               % "Length"
                                               % "InputConcent"
                                               % "Observed").str());
        for (const auto &i : ls)
        {
            Locus l;
            
            switch (o.metrs)
            {
                case RExpress::Metrics::Gene:    { l = r.findGene(ChrIS, i.first)->l;  break; }
                case RExpress::Metrics::Isoform: { l = r.findTrans(ChrIS, i.first)->l; break; }
            }

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
        o.info("Genome: " + toString(i.gData.size()));
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
     * Generating RnaExpression_queries.csv
     */
    
    writeQueries("RnaExpression_queries.csv", stats, o);
    
    /*
     * Generating RnaExpression_express.R
     */
    
    RExpress::generateR("RnaExpression_express.R", "RnaExpression_sequins.csv", stats, o);
}