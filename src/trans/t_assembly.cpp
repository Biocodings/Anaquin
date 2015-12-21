#include <fstream>
#include "data/compare.hpp"
#include "trans/t_assembly.hpp"
#include "parsers/parser_gtf.hpp"

#define CHECK_AND_SORT(t) { assert(!t.empty()); std::sort(t.begin(), t.end(), [](const Feature& x, const Feature& y) { return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end); }); }

using namespace Anaquin;

// Defined for cuffcompare
Compare __cmp__;

// Defined for cuffcompare
extern int cuffcompare_main(const char *ref, const char *query);

template <typename F> static void extractIntrons(const std::map<SequinID, std::vector<Feature>> &x, F f)
{
    Feature ir;

    for (const auto & ts : x)
    {
        for (auto i = 0; i < ts.second.size(); i++)
        {
            if (i)
            {
                if (ts.second[i-1].tID == ts.second[i].tID)
                {
                    ir = ts.second[i];
                    ir.l = Locus(ts.second[i - 1].l.end + 1, ts.second[i].l.start - 1);
                    f(ts.second[i-1], ts.second[i], ir);
                }
            }
        }
    }
}

static std::string createFilteredGTF(const FileName &file, const ChromoID &cID)
{
    std::string line;
    const auto tmp = tmpnam(NULL);
    
    std::ofstream out(tmp);

    ParserGTF::parse(file, [&](const Feature &f, const std::string &l, const ParserProgress &)
    {
        if (f.id == cID)
        {
            out << l << std::endl;
        }
    });
    
    out.close();
    
    return tmp;
}

static std::string summary()
{
    return "Summary for dataset: %1%\n\n"
           "   Experiment: %2% features\n"
           "   Synthetic:  %3% features\n\n"
           "   Reference:  %4% exons\n"
           "   Reference:  %5% introns\n\n"
           "   ***\n"
           "   *** The following statistics are computed for exact and fuzzy.\n"
           "   ***\n"
           "   *** The fuzzy level is 10 nucleotides.\n"
           "   ***\n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %6% (%7%)\n"
           "   Specificity: %8% (%9%)\n"
           "   Detection:   %10% (%11%)\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %12% (%13%)\n"
           "   Specificity: %14% (%15%)\n"
           "   Detection:   %16% (%17%)\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %18%\n"
           "   Specificity: %19%\n"
           "   Detection:   %20% (%21%)\n\n"
           "   -------------------- Intron Chain level --------------------\n\n"
           "   Sensitivity: %38% (%39%)\n"
           "   Specificity: %40% (%41%)\n\n"
           "   -------------------- Transcript level --------------------\n\n"
           "   Sensitivity: %22% (%23%)\n"
           "   Specificity: %24% (%25%)\n\n"
           "   Missing exons: %26%/%27% (%28%)\n"
           "   Missing introns: %29%/%30% (%31%)\n\n"
           "   Novel exons: %32%/%33% (%34%)\n"
           "   Novel introns: %35%/%36% (%37%)\n\n";
}

TAssembly::Stats TAssembly::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_trans;

    assert(!o.ref.empty() && !o.query.empty());
    
    /*
     * 1. Initalize the statistics
     */
    
    TAssembly::Stats stats;
    
    stats.data[ChrT];
    
    if (r.chromoIDs().size() > 1)
    {
        stats.data[ChrE];
    }
    
    stats.eHist = r.hist();
    stats.iHist = r.hist();
    stats.tHist = r.hist();
    stats.bHist = r.geneHist(ChrT);

    /*
     * 2. Filtering transcripts
     */
    
    auto copyStats = [&](const ChromoID &cID)
    {
        stats.data[cID].eSN  = std::min(__cmp__.e_sn  / 100.0, 1.0);
        stats.data[cID].eSP  = std::min(__cmp__.e_sp  / 100.0, 1.0);
        stats.data[cID].eFSN = std::min(__cmp__.e_fsn / 100.0, 1.0);
        stats.data[cID].eFSP = std::min(__cmp__.e_fsp / 100.0, 1.0);

        stats.data[cID].iSN  = std::min(__cmp__.i_sn  / 100.0, 1.0);
        stats.data[cID].iSP  = std::min(__cmp__.i_sp  / 100.0, 1.0);
        stats.data[cID].iFSN = std::min(__cmp__.i_fsn / 100.0, 1.0);
        stats.data[cID].iFSP = std::min(__cmp__.i_fsp / 100.0, 1.0);
        
        stats.data[cID].tSN  = std::min(__cmp__.t_sn  / 100.0, 1.0);
        stats.data[cID].tSP  = std::min(__cmp__.t_sp  / 100.0, 1.0);
        stats.data[cID].tFSN = std::min(__cmp__.t_fsn / 100.0, 1.0);
        stats.data[cID].tFSP = std::min(__cmp__.t_fsp / 100.0, 1.0);
        
        stats.data[cID].bSN  = std::min(__cmp__.b_sn  / 100.0, 1.0);
        stats.data[cID].bSP  = std::min(__cmp__.b_sp  / 100.0, 1.0);
    };
    
    auto compareGTF = [&](const ChromoID &cID)
    {
        const auto query = createFilteredGTF(file, cID);
        o.logInfo("Filtered transcript: " + query + " has been created");
        
        o.logInfo("Invoking Cuffcompare: " + o.ref);
        o.logInfo("Invoking Cuffcompare: " + query);
        
        if (cuffcompare_main(o.ref.c_str(), query.c_str()))
        {
            throw std::runtime_error("Failed to analyze the given transcript. Please check the file and try again.");
        }
    };

    o.info("Generating a filtered transcript");

    std::for_each(stats.data.begin(), stats.data.end(), [&](const std::pair<ChromoID, TAssembly::Stats::Data> &p)
    {
        compareGTF(p.first);
        copyStats(p.first);
    });
    
    /*
     * 3. Classifying the transcript (only for chrT because detection limit is needed)
     */

    o.info("Parsing transcript");
    
    std::vector<Feature> q_exons;
    std::map<SequinID, std::vector<Feature>> q_exons_;

    Confusion t;
    
    ParserGTF::parse(file, [&](const Feature &f, const std::string &, const ParserProgress &p)
    {
        if (!(p.i % 1000000)) { o.wait(std::to_string(p.i)); }
        
        if (f.id != Standard::chrT)
        {
            stats.n_expT++;
        }
        else
        {
            stats.n_chrT++;
        }

        auto classifyT = [&]()
        {
            const auto &r = Standard::instance().r_trans;
            
            if (f.cID != ChrT)
            {
                return;
            }
            
            switch (f.type)
            {
                case Exon:
                {
                    const TransRef::ExonData *match;
                    
                    q_exons.push_back(f);
                    q_exons_[f.tID].push_back(f);
                    
                    if (classify(t, f, [&](const Feature &)
                    {
                        return (match = r.findExon(ChrT, f.l, Exact));
                    }))
                    {
                        stats.eHist.at(match->iID)++;
                    }
                    
                    break;
                }
                    
                case Transcript:
                {
                    const TransData *match;
                    
                    if (classify(t, f, [&](const Feature &)
                    {
                        return (match = r.match(f.l, Overlap));
                    }))
                    {
                        stats.tHist.at(match->id)++;
                    }
                    
                    break;
                }
                    
                // There're many other possibilties in a GTF file, but we don't need them
                default: { break; }
            }
        };

        classifyT();
    });
    
    /*
     * 3. Generating introns
     */
    
    o.info("Generating introns");
    
    /*
     * Sort the query exons as there is no guarantee that those are sorted
     */
    
    for (auto &i : q_exons_)
    {
        CHECK_AND_SORT(i.second);
    }
    
    /*
     * Now that the query exons are sorted. We can extract and classify the introns for each pair
     * of successive exon.
     */
    
    const TransRef::IntronData *match;
    
    extractIntrons(q_exons_, [&](const Feature &, const Feature &, Feature &i)
    {
        if (classify(t, i, [&](const Feature &)
        {
            return (match  = r.findIntron("chrT", i.l, Exact));
        }))
        {
            stats.iHist.at(match->iID)++;
        }
    });
    
    /*
     * The counts for query bases is the total non-overlapping bases of all the exons in the experiment.
     */
    
    //countBase(r.mergedExons(), q_exons, t, stats.chrT->hb);
    
    /*
     * 4. Collecting statistics
     */
    
    o.info("Collecting statistics");
    
    stats.eLimit = r.limit(stats.eHist);
    stats.tLimit = r.limit(stats.tHist);
    stats.iLimit = r.limit(stats.iHist);
    stats.bLimit = r.limitGene(stats.bHist);
    
    return stats;
}

static void writeSummary(const FileName &file,
                         const FileName &src,
                         const TAssembly::Stats &stats,
                         const TAssembly::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    o.writer->open(file);
    o.writer->write((boost::format(summary()) % file
                                              % stats.n_expT
                                              % stats.n_chrT
                                              % r.data().size()
                                              % r.countIntrons(ChrT)
                                              % (__cmp__.e_sn  / 100.0) // 6
                                              % (__cmp__.e_fsn / 100.0)
                                              % (__cmp__.e_sp  / 100.0)
                                              % (__cmp__.e_fsp / 100.0)
                                              % (stats.eLimit.id.empty() ? "-" : std::to_string(stats.eLimit.abund))
                                              %  stats.eLimit.id
                                              % (__cmp__.i_sn  / 100.0) // 12
                                              % (__cmp__.i_fsn / 100.0)
                                              % (__cmp__.i_sp  / 100.0)
                                              % (__cmp__.i_fsp / 100.0)
                                              % (stats.iLimit.id.empty() ? "-" : std::to_string(stats.iLimit.abund))
                                              %  stats.iLimit.id
                                              % (__cmp__.b_sn / 100.0) // 18
                                              % (__cmp__.b_sp / 100.0)
                                              % (stats.bLimit.id.empty() ? "-" : std::to_string(stats.bLimit.abund)) // 20
                                              %  stats.bLimit.id
                                              % (__cmp__.t_sn  / 100.0)
                                              % (__cmp__.t_fsn / 100.0)
                                              % (__cmp__.t_sp  / 100.0)
                                              % (__cmp__.t_fsp / 100.0) // 25
                                              % (__cmp__.missedExonsN)
                                              % (__cmp__.missedExonsR)
                                              % (__cmp__.missedExonsP / 100.0)
                                              % (__cmp__.missedIntronsN)
                                              % (__cmp__.missedIntronsR)         // 30
                                              % (__cmp__.missedIntronsP / 100.0)
                                              % (__cmp__.novelExonsN)
                                              % (__cmp__.novelExonsR)
                                              % (__cmp__.novelExonsP / 100.0)
                                              % (__cmp__.novelIntronsN)
                                              % (__cmp__.novelIntronsR)          // 36
                                              % (__cmp__.novelIntronsP / 100.0)
                                              % (__cmp__.c_sn  / 100.0)
                                              % (__cmp__.c_fsn / 100.0)
                                              % (__cmp__.c_sp  / 100.0)
                                              % (__cmp__.c_fsp / 100.0)).str());
    o.writer->close();
}

static void writeSequins(const FileName &file, const FileName &src, const TAssembly::Stats &stats, const TAssembly::Options &o)
{
    o.writer->open(file);
    o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());
    
    auto format = "%1%\t%2%\t%3%\t%4%";
    o.writer->write((boost::format(format) % "ID" % "Exon" % "Intron" % "Transcript").str());
    
    for (const auto &i : stats.eHist)
    {
        o.writer->write((boost::format(format) % i.first
                                               % stats.eHist.at(i.first)
                                               % stats.iHist.at(i.first)
                                               % stats.tHist.at(i.first)).str());
    }

    /*
     o.writer->write("\n");
     
     format = "%1%\t%2%";
     o.writer->write((boost::format(format) % "ID" % "Base").str());
     
     for (const auto &i : stats.hb)
     {
     o.writer->write((boost::format(format) % i.first % stats.hb.at(i.first)).str());
     }
     */

    o.writer->close();
}

void TAssembly::report(const FileName &file, const Options &o)
{
    const auto stats = TAssembly::analyze(file, o);

    o.info("Generating summary statistics");
    writeSummary("TransAssembly_summary.stats", file, stats, o);

    o.info("Generating sequin statistics");
    writeSequins("TransAssembly_quins.stats", file, stats, o);
}

void TAssembly::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TAssembly::analyze(files, o);

    /*
     * Generating summary statistics for each replicate
     */
    
    o.info("Generating summary statistics");

    for (auto i = 0; i < files.size(); i++)
    {
        writeSummary((boost::format("TransAssembly_%1%_summary.stats") % files[i]).str(), files[i], stats[i], o);
    }
    
    /*
     * Generating sequin statistics for each replicate
     */
    
    for (auto i = 0; i < files.size(); i++)
    {
        writeSequins((boost::format("TransAssembly_%1%_sequin.stats") % files[i]).str(), files[i], stats[i], o);
    }
}