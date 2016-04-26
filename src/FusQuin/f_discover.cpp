#include "FusQuin/FUSQuin.hpp"
#include "FusQuin/f_discover.hpp"

using namespace Anaquin;

extern Scripts PlotFROC();

template <typename T> ChrID toChrTEndo(const T &t)
{
    return t == ChrT ? ChrT : Geno;
}

FDiscover::Stats FDiscover::analyze(const FileName &file, const FDiscover::Options &o)
{
    const auto &r = Standard::instance().r_fus;
    
    FDiscover::Stats stats;
    
    stats.data[ChrT];
    stats.data[ChrT].hist = r.fusionHist();
    stats.data[Geno];

    FUSQuin::analyze<FDiscover::Options>(file, o, [&](const FUSQuin::Match &match)
    {
        if (match.label == FUSQuin::Label::Positive)
        {
            const auto cID = toChrTEndo(match.query.cID_1);

            stats.data[cID].tps.push_back(match);
            stats.data[cID].hist[match.known->id]++;
        }
        else
        {
            if (match.query.cID_1 == ChrT && match.query.cID_2 == ChrT)
            {
                stats.data[ChrT].fps.push_back(match);
            }
            else
            {
                stats.data[Geno].fps.push_back(match);
            }
        }
    });
    
    /*
     * Find out all the missing references
     */

    for (auto &i : stats.data)
    {
        for (const auto &j : i.second.hist)
        {
            if (i.first == ChrT)
            {
                i.second.fns.push_back(*r.findFusion(j.first));
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
                         "   *** Number of fusions detected in the synthetic chromosome\n"
                         "   ***\n\n"
                         "   Synthetic: %2% fusions\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %3%\n\n"
                         "   Synthetic: %4% fusions\n\n"
                         "   ************************************************************\n"
                         "   ***                                                      ***\n"
                         "   ***        Statistics for the synthetic chromosome       ***\n"
                         "   ***                                                      ***\n"
                         "   ************************************************************\n\n"
                         "   True Positives:  %5% fusions\n"
                         "   False Positives: %6% fusions\n\n"
                         "   ***\n"
                         "   *** Performance metrics\n"
                         "   ***\n\n"
                         "   Sensitivity: %7%\n"
                         "   Specificity: %8%\n\n";
    
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % stats.countDetect(ChrT)
                                            % o.rChrT
                                            % r.countFusion()
                                            % stats.countTP(ChrT)
                                            % stats.countFP(ChrT)
                                            % stats.sn(ChrT)
                                            % stats.pc(ChrT)).str());
    o.writer->close();
}

static void writeClass(const FileName &file, const ChrID &cID, const FDiscover::Stats &stats, const FDiscover::Options &o)
{
    const auto &data  = stats.data.at(cID);
    const auto format = "%1%\t%2%\t%3%\t%4%";

    o.writer->open(file);
    o.writer->write((boost::format(format) % "Sequin"
                                           % "Label"
                                           % "Position_1"
                                           % "Position_2").str());

    for (const auto &tp : data.tps)
    {
        o.writer->write((boost::format(format) % tp.known->id
                                               % "TP"
                                               % tp.query.l1
                                               % tp.query.l2).str());
    }

    for (const auto &fp : data.fps)
    {
        o.writer->write((boost::format(format) % "-"
                                               % "FP"
                                               % fp.query.l1
                                               % fp.query.l2).str());
    }
}

static void writeQuins(const FileName &file, const FDiscover::Stats &stats, const FDiscover::Options &o)
{
    o.writer->open(file);
    
    const auto format = "%1%\t%2%";
    
    o.writer->write((boost::format(format) % "ID"
                                           % "Counts").str());

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
     * Generating summary statistics
     */

    writeSummary("FusionDiscover_summary.stats", stats, o);

    /*
     * Generating classified statistics for the fusions
     */
    
    writeClass("FusionDiscover_labels.csv", ChrT, stats, o);
    
    /*
     * Generating ROC curve
     */
    
    o.writer->open("FusionDiscover_ROC.R");
    o.writer->write(RWriter::createScript("FusionDiscover_labels.csv", PlotFROC()));
    o.writer->close();

    /*
     * Generating sequin statistics
     */

    writeQuins("FusionDiscover_quins.stats", stats, o);
}