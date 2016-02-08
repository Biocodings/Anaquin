#include "fusion/FUSQuin.hpp"
#include "fusion/f_discover.hpp"

using namespace Anaquin;

FDiscover::Stats FDiscover::analyze(const FileName &file, const FDiscover::Options &o)
{
    FDiscover::Stats stats;
    
    stats.data[ChrT];
    stats.data[ChrT].hist = Standard::instance().r_fus.hist();
    stats.data[Endo];

    FUSQuin::analyze<FDiscover::Options>(file, o, [&](const FUSQuin::Match &match)
    {
        stats.hist[match.query.cID_1]++;
        stats.hist[match.query.cID_2]++;

        if (match.state == FUSQuin::State::Positive)
        {
            if (match.query.cID_1 == ChrT)
            {
                stats.data[ChrT].tps.push_back(match.known);
            }
            else
            {
                stats.data[Endo].tps.push_back(match.known);
            }
        }
        else
        {
            if (match.query.cID_1 == ChrT && match.query.cID_2 == ChrT)
            {
                stats.data[ChrT].fps.push_back(match.query);
            }
            else
            {
                stats.data[Endo].fps.push_back(match.query);
            }
        }
    });
    
    return stats;
}

static void writeSummary(const FileName &file, const FDiscover::Stats &stats, const FDiscover::Options &o)
{
    o.writer->open("FusionDiscover_summary.stats");
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   Experiment: %2% fusions\n"
                         "   Synthetic: %3% fusions\n"
                         "   Synthetic-Experiment: %4% fusions\n"
                         "   Reference: %5% sequins\n\n"
                         "   Fuzzy: %6%\n\n"
                         "   Sensitivity: %7%\n"
                         "   Specificity: %8%\n";
/*
    o.writer->write((boost::format(summary) % file
                     % stats.chrT->n_endo
                     % stats.chrT->n_chrT
                     % stats.chrT->hg38_chrT
                     % stats.chrT->m.nr()
                     % o.fuzzy
                     % stats.chrT->m.sn()
                     % stats.chrT->m.ac()).str());
*/

    o.writer->close();    
}

void FDiscover::report(const FileName &file, const FDiscover::Options &o)
{
    const auto stats = analyze(file, o);

    /*
     * Generating summary statistics
     */
    
    o.info("Generating summary statistics");

    /*
     * Generating labels for the known fusions
     */
    
    /*
     * Generating detailed statistics for sequins
     */

//
//
//    {
//        o.info("Generating sequins statistics");
//        o.writer->open("FusionDiscover_quins.stats");
//
//        const auto format = "%1%\t%2%";
//
//        o.writer->write((boost::format(format) % "ID" % "Counts").str());
//
//        for (const auto &i : stats.chrT->h)
//        {
//                o.writer->write((boost::format(format) % i.first
//                                                       % i.second).str());
//        }
//
//        o.writer->close();
//    }
}