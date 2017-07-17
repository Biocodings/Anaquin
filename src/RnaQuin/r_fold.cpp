/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
 */

#include "tools/tools.hpp"
#include "RnaQuin/r_fold.hpp"
#include "tools/gtf_data.hpp"
#include "RnaQuin/RnaQuin.hpp"
#include "parsers/parser_fold.hpp"
#include "parsers/parser_edgeR.hpp"
#include "parsers/parser_cdiff.hpp"
#include "parsers/parser_sleuth.hpp"
#include "parsers/parser_DESeq2.hpp"

using namespace Anaquin;

extern Scripts PlotTROC();

typedef RFold::Metrics Metrics;

static bool shouldAggregate(const RFold::Options &o)
{
    return o.metrs == Metrics::Gene && o.format == RFold::Format::Sleuth;
}

std::vector<std::string> RFold::classify(const std::vector<double> &qs, const std::vector<double> &folds, double qCut, double foldCut)
{
    assert(qs.size() == folds.size());
    
    std::vector<std::string> r;
    
    for (auto i = 0; i < qs.size(); i++)
    {
        // Differential expressed
        if (qs[i] <= qCut)
        {
            r.push_back(fabs(folds[i]) <= foldCut ? "FP" : "TP");
        }
        
        // Non-differential expressed
        else
        {
            r.push_back(fabs(folds[i]) <= foldCut ? "TN" : "FN");
        }
    }
    
    return r;
}

template <typename T> void classifySyn(RFold::Stats &stats, const T &t, Metrics metrs, const RFold::Options &o)
{
    if (t.status == DiffTest::Status::NotTested)
    {
        return;
    }
    
    const auto &r = Standard::instance().r_rna;

    auto f = [&](const SequinID &id, Concent exp)
    {
        if (t.p == 0) { o.warn(id + " gives p-value of 0"); }
        if (t.q == 0) { o.warn(id + " gives q-value of 0"); }

        stats.data[id].p     = t.p;
        stats.data[id].q     = t.q;
        stats.data[id].exp   = exp;
        stats.data[id].obs   = t.logF_;
        stats.data[id].se    = t.logFSE;
        stats.data[id].mean  = t.mean;
        stats.data[id].samp1 = t.samp1;
        stats.data[id].samp2 = t.samp2;

        if (!isnan(exp) && !isnan(t.logF_) && std::isfinite(t.logF_))
        {
            stats.add(id, exp, t.logF_);
            stats.nSeqs++;
        }
    };
    
    switch (metrs)
    {
        case Metrics::Gene:
        {
//            assert(!t.gID.empty());
//      TODO      const auto match = r.findGene(t.cID, t.gID);
//            
//            if (match)
//            {
//                stats.hist.at(t.cID).at(t.gID)++;
//                f(t.gID, r.logFoldGene(t.gID));
//            }
//            else
//            {
//                o.warn(t.gID + " not found");
//            }

            break;
        }

        case Metrics::Isoform:
        {
            assert(!t.iID.empty());
            const auto match = r.input1(t.iID); // r.match(t.iID);
            
            if (match)
            {
                // dont need this line stats.hist.at(t.cID).at(t.iID)++;
                // TOODf(t.iID, r.logFoldSeq(t.iID));
            }
            else
            {
                o.warn(t.iID + " not found");
            }
            
            break;
        }
    }
}

template <typename T> void update(RFold::Stats &stats, const T &x, Metrics metrs, const RFold::Options &o)
{
    typedef DiffTest::Status Status;
    
    if (isRNARevChr(x.cID))
    {
        classifySyn(stats, x, metrs, o);
    }
    else
    {
        stats.nEndo++;
    }
}

template <typename Functor> RFold::Stats calculate(const RFold::Options &o, Functor f)
{
    const auto &r = Standard::instance().r_rna;

    RFold::Stats stats;

    switch (o.metrs)
    {
        case Metrics::Gene:    { stats.hist = r.histGene(); break; }
        case Metrics::Isoform: { stats.hist = r.histIsof(); break; }
    }

    assert(!stats.hist.empty());

    f(stats);
    
    return stats;
}

RFold::Stats RFold::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    return calculate(o, [&](RFold::Stats &stats)
    {
        switch (o.format)
        {
            case Format::Sleuth:
            {
                ParserSleuth::parse(file, [&](const ParserSleuth::Data &x, const ParserProgress &)
                {
                    // Should we aggregate because this is at the gene level?
                    if (shouldAggregate(o))
                    {
                        update(stats, x, RFold::Metrics::Isoform, o);
                    }
                    else
                    {
                        update(stats, x, o.metrs, o);
                    }
                });

                break;
            }
                
            case Format::Anaquin:
            {
                ParserDiff::parse(file, [&](const ParserDiff::Data &x, const ParserProgress &)
                {
                    update(stats, x, o.metrs, o);
                });

                break;
            }

            case Format::DESeq2:
            {
                ParserDESeq2::parse(file, [&](const ParserDESeq2::Data &x, const ParserProgress &)
                {
                    update(stats, x, o.metrs, o);
                });
                
                break;
            }
                
            case Format::edgeR:
            {
                ParserEdgeR::parse(file, [&](const ParserEdgeR::Data &x, const ParserProgress &)
                {
                    update(stats, x, o.metrs, o);
                });
                
                break;
            }
                
            case Format::Cuffdiff:
            {
                ParserCDiff::parse(file, [&](const ParserCDiff::Data &x, const ParserProgress &)
                {
                    update(stats, x, o.metrs, o);
                });

                break;
            }
        }
        
        if (shouldAggregate(o))
        {
            const auto &r = Standard::instance().r_rna;

            std::map<GeneID, Fold> g2f;
            
            for (const auto &i : stats.data)
            {
                const auto m = r.findTrans(ChrIS, i.first);
                
                // We should always find the gene identifier
                A_ASSERT(m);
                
                g2f[m->gID] += i.second.obs;
            }
            
            stats = RFold::Stats();
            
            for (const auto &i : g2f)
            {
                stats.nSeqs++;
                stats.data[i.first].obs = i.second;
                stats.data[i.first].exp = r.input2(i.first);
                stats.add(i.first, i.second, r.input2(i.first));
            }
        }
    });
}

void RFold::writeCSV(const FileName &file, const RFold::Stats &stats, const RFold::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RFold::generateCSV(stats, o));
    o.writer->close();
}

Scripts RFold::generateCSV(const RFold::Stats &stats, const RFold::Options &o)
{
    const auto &r = Standard::instance().r_rna;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%";

    std::stringstream ss;
    ss << (boost::format(format) % "ID"
                                 % "Length"
                                 % "Sample1"
                                 % "Sample2"
                                 % "ExpLFC"
                                 % "ObsLFC"
                                 % "SD"
                                 % "Pval"
                                 % "Qval"
                                 % "Mean").str();
    
    const auto ids = o.metrs == RFold::Metrics::Gene ? r.seqsL2() : r.seqsL1();

    // For each sequin gene or isoform...
    for (const auto id : ids)
    {
        Locus l;
        LogFold fold;

        switch (o.metrs)
        {
            case Metrics::Gene:
            {
                fold = r.input2(id);
                l = r.findGene(ChrIS, id)->l;
                break;
            }

            case Metrics::Isoform:
            {
                fold = r.input1(id);
                l = r.findTrans(ChrIS, id)->l;
                break;
            }
        }

        // Undetected or nothing informative ...
        if (!stats.data.count(id) || isnan(stats.data.at(id).obs))
        {
            ss << (boost::format(format) % id
                                         % l.length()
                                         % "NA"
                                         % "NA"
                                         % toString(fold)
                                         % "NA"
                                         % "NA"
                                         % "NA"
                                         % "NA"
                                         % "NA").str() << std::endl;
            continue;
        }
        
        const auto &x = stats.data.at(id);
        A_ASSERT(fold == x.exp);
        
        ss << (boost::format(format) % id
                                     % l.length()
                                     % toString(x.samp1)
                                     % toString(x.samp2)
                                     % toString(fold)
                                     % toString(x.obs)
                                     % toString(x.se)
                                     % ld2ss(x.p)
                                     % ld2ss(x.q)
                                     % toString(x.mean)).str() << std::endl;
    }

    return ss.str();
}

void RFold::writeSummary(const FileName &file,
                            const FileName &src,
                            const RFold::Stats &stats,
                            const RFold::Options &o,
                            const Units &units)
{
    o.writer->open(file);
    o.writer->write(RFold::generateSummary(src, stats, o, units));
    o.writer->close();
}

Scripts RFold::generateSummary(const FileName &src,
                               const RFold::Stats &stats,
                               const RFold::Options &o,
                               const Units &units)
{
    const auto &r = Standard::instance().r_rna;
    const auto lm = stats.linear(false);
    
    // No reference coordinate annotation given here
    const auto nSyn = o.metrs == Metrics::Gene ? r.seqsL2().size() : r.seqsL1().size();
    
    const auto title = (o.metrs == Metrics::Gene ? "Genes Expressed" : "Isoform Expressed");
    
    const auto summary = "-------RnaFoldChange Output\n\n"
                         "       Summary for input: %1%\n\n"
                         "-------Reference Annotations\n\n"
                         "       Synthetic: %2% %3%\n"
                         "       Mixture file: %4%\n\n"
                         "-------%5%\n\n"
                         "       Synthetic: %6% %3%\n"
                         "       Genome:    %7% %3%\n\n"
                         "-------Linear regression (log2 scale)\n\n"
                         "       Slope:       %8%\n"
                         "       Correlation: %9%\n"
                         "       R2:          %10%\n"
                         "       F-statistic: %11%\n"
                         "       P-value:     %12%\n"
                         "       SSM:         %13%, DF: %14%\n"
                         "       SSE:         %15%, DF: %16%\n"
                         "       SST:         %17%, DF: %18%\n";

    return (boost::format(summary) % src         // 1
                                   % nSyn        // 2
                                   % units       // 3
                                   % LadRef()    // 4
                                   % title       // 5
                                   % stats.nSeqs // 6
                                   % stats.nEndo // 7
                                   % lm.m        // 8
                                   % lm.r        // 9
                                   % lm.R2       // 10
                                   % lm.F        // 11
                                   % lm.p        // 12
                                   % lm.SSM      // 13
                                   % lm.SSM_D    // 14
                                   % lm.SSE      // 15
                                   % lm.SSE_D    // 16
                                   % lm.SST      // 17
                                   % lm.SST_D    // 18
                     ).str();
}

Scripts RFold::generateRFold(const RFold::Stats &stats, const FileName &csv, const RFold::Options &o)
{
    //const auto extra = o.format == Format::DESeq2 ? ", std=data$SD" : "";
    const auto extra = o.format == Format::DESeq2 ? "" : "";

    return RWriter::createFold(csv,
                               o.work,
                               o.metrs == RFold::Metrics::Gene ? "Gene Fold Change" : "Isoform Fold Change",
                               "Expected fold change (log2)",
                               "Measured fold change (log2)",
                               "ExpLFC",
                               "ObsLFC",
                               false,
                               extra);
}

void RFold::writeRFold(const FileName &file, const RFold::Stats &stats, const RFold::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RFold::generateRFold(stats, "RnaFoldChange_sequins.csv", o));
    o.writer->close();
}

Scripts RFold::generateRROC(const RFold::Stats &stats, const RFold::Options &o)
{
    return RWriter::createScript("RnaFoldChange_sequins.csv", PlotTROC());
}

void RFold::writeRROC(const FileName &file, const RFold::Stats &stats, const RFold::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RFold::generateRROC(stats, o));
    o.writer->close();
}

Scripts RFold::generateRLODR(const RFold::Stats &stats, const RFold::Options &o)
{
    return RWriter::createScript("RnaFoldChange_sequins.csv", PlotTLODR());
}

void RFold::writeRLODR(const FileName &file, const RFold::Stats &stats, const RFold::Options &o)
{
    switch (o.format)
    {
        case Format::edgeR:
        case Format::Cuffdiff:
        {
            o.info("Skip RnaFoldChange_LODR.R because no average counts given");
            break;
        }

        case Format::Sleuth:
        case Format::DESeq2:
        case Format::Anaquin:
        {
            o.generate(file);
            o.writer->open(file);
            o.writer->write(RFold::generateRLODR(stats, o));
            o.writer->close();
            break;
        }
    }
}

void RFold::report(const FileName &file, const Options &o)
{
    const auto m = std::map<Metrics, std::string>
    {
        { Metrics::Gene,    "genes"    },
        { Metrics::Isoform, "isoforms" },
    };

    switch (o.metrs)
    {
        case Metrics::Gene:    { o.info("Gene Differential");    break; }
        case Metrics::Isoform: { o.info("Isoform Differential"); break; }
    }
    
    const auto stats = RFold::analyze(file, o);
    const auto units = m.at(o.metrs);
    
    o.info("Generating statistics");
    
    /*
     * Generating RnaFoldChange_summary.stats
     */

    RFold::writeSummary("RnaFoldChange_summary.stats", file, stats, o, units);

    /*
     * Generating RnaFoldChange_sequins.csv
     */

    RFold::writeCSV("RnaFoldChange_sequins.csv", stats, o);
    
    /*
     * Generating RnaFoldChange_fold.R
     */
    
    RFold::writeRFold("RnaFoldChange_fold.R", stats, o);
    
    /*
     * Generating RnaFoldChange_ROC.R
     */

    RFold::writeRROC("RnaFoldChange_ROC.R", stats, o);

    /*
     * Generating RnaFoldChange_LODR.R
     */
    
    RFold::writeRLODR("RnaFoldChange_LODR.R", stats, o);
}
