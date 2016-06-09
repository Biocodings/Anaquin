#include "FusQuin/FUSQuin.hpp"
#include "FusQuin/f_discover.hpp"

using namespace Anaquin;

extern Scripts PlotFROC();

FDiscover::Stats FDiscover::analyze(const FileName &file, const FDiscover::Options &o)
{
    const auto &r = Standard::instance().r_fus;
    
    FDiscover::Stats stats;
    
    stats.data[Geno];
    stats.data[ChrT].hist = r.fusionHist();

    FUSQuin::analyze<FDiscover::Options>(file, o, [&](const FUSQuin::Match &match)
    {
        switch (match.label)
        {
            case FUSQuin::Label::Positive:
            {
                stats.data[ChrT].tps.push_back(match);
                stats.data[ChrT].hist[match.known->id]++;
                break;
            }
                
            case FUSQuin::Label::GenoChrT:
            case FUSQuin::Label::Negative:
            {
                stats.data[ChrT].fps.push_back(match);
                break;
            }
                
            case FUSQuin::Label::Geno:
            {
                stats.data[Geno].fps.push_back(match);
                break;
            }
        }
    });
    
    /*
     * Find out all the missing fusions (only chrT for now)
     */

    for (auto &i : stats.data)
    {
        if (i.first == ChrT)
        {
            for (const auto &j : i.second.hist)
            {
                if (!j.second)
                {
                    i.second.fns.push_back(*r.findFusion(j.first));
                }
            }
        }
    }

    return stats;
}

static void writeSummary(const FileName &file, const FDiscover::Stats &stats, const FDiscover::Options &o)
{
    const auto &r = Standard::instance().r_fus;

    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Number of fusions detected in the synthetic and genome\n"
                         "   ***\n\n"
                         "   Synthetic: %2% fusions\n"
                         "   Genome:    %3% fusions\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %4%\n\n"
                         "   Synthetic: %5% fusions\n\n"
                         "   Fuzzy:     %6%\n\n"
                         "   ***                                         \n"
                         "   *** Statistics for the synthetic chromosome \n"
                         "   ***                                         \n\n"
                         "   True Positives:  %7% fusions\n"
                         "   False Positives: %8% fusions\n\n"
                         "   ***\n"
                         "   *** Performance metrics\n"
                         "   ***\n\n"
                         "   Sensitivity: %9%\n"
                         "   Precision:   %10%\n\n";
    
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % stats.countDetect(ChrT)
                                            % stats.countDetect(Geno)
                                            % o.rAnnot
                                            % r.countFusion()
                                            % o.fuzzy
                                            % stats.countTP(ChrT)
                                            % stats.countFP(ChrT)
                                            % stats.sn(ChrT)
                                            % stats.pc(ChrT)).str());
    o.writer->close();
}

static void writeQuery(const FileName &file, const ChrID &cID, const FDiscover::Stats &stats, const FDiscover::Options &o)
{
    const auto &data  = stats.data.at(cID);
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";

    o.writer->open(file);
    o.writer->write((boost::format(format) % "seq"
                                           % "label"
                                           % "chr1"
                                           % "chr2"
                                           % "str1"
                                           % "str2"
                                           % "pos1"
                                           % "pos2"
                                           % "reads").str());

    for (const auto &tp : data.tps)
    {
        o.writer->write((boost::format(format) % tp.known->id
                                               % "TP"
                                               % tp.query.cID_1
                                               % tp.query.cID_2
                                               % tp.query.s1
                                               % tp.query.s2
                                               % tp.query.l1
                                               % tp.query.l2
                                               % tp.query.reads).str());
    }

    for (const auto &fp : data.fps)
    {
        o.writer->write((boost::format(format) % "-"
                                               % "FP"
                                               % fp.query.cID_1
                                               % fp.query.cID_2
                                               % fp.query.s1
                                               % fp.query.s2
                                               % fp.query.l1
                                               % fp.query.l2
                                               % fp.query.reads).str());
    }
}

static void writeQuins(const FileName &file, const FDiscover::Stats &stats, const FDiscover::Options &o)
{
    o.writer->open(file);
    
    const auto format = "%1%\t%2%";
    
    o.writer->write((boost::format(format) % "seq"
                                           % "reads").str());

    for (const auto &i : stats.data.at(ChrT).hist)
    {
        o.writer->write((boost::format(format) % i.first
                                               % i.second).str());
    }
    
    o.writer->close();
}

void FDiscover::report(const FileName &file, const FDiscover::Options &o)
{
    const auto stats = analyze(file, o);

    o.info("Generating statistics");

    /*
     * Generating FusDiscover_summary.stats
     */

    writeSummary("FusDiscover_summary.stats", stats, o);

    /*
     * Generating FusDiscover_quins.stats
     */
    
    writeQuins("FusDiscover_quins.stats", stats, o);
    
    /*
     * Generating FusDiscover_queries.stats
     */
    
    writeQuery("FusDiscover_queries.stats", ChrT, stats, o);
    
    /*
     * Generating FusDiscover_ROC.R
     */
    
    o.writer->open("FusDiscover_ROC.R");
    o.writer->write(RWriter::createScript("FusDiscover_query.stats", PlotFROC()));
    o.writer->close();
    o.report->addFile("FusDiscover_ROC.R", "FusDiscover_ROC.R");
    
    /*
     * Generating FusDiscover_report.pdf
     */

    o.report->open("FusDiscover_report.pdf");
    o.report->addTitle("FusDiscover_report");
    o.report->addFile("FusDiscover_summary.stats");
    o.report->addFile("FusDiscover_quins.stats");
    o.report->addFile("FusDiscover_query.stats");
    o.report->addFile("FusDiscover_ROC.R");
}