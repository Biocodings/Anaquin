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
           "   Sensitivity: %22% (%23%)\n"
           "   Specificity: %24% (%25%)\n\n"
           "   -------------------- Transcript level --------------------\n\n"
           "   Sensitivity: %26% (%27%)\n"
           "   Specificity: %28% (%29%)\n\n"
           "   Missing exons: %30%/%31% (%32%)\n"
           "   Missing introns: %33%/%34% (%35%)\n\n"
           "   Novel exons: %36%/%37% (%38%)\n"
           "   Novel introns: %39%/%40% (%41%)\n\n";
}

static TAssembly::Stats init()
{
    const auto &r = Standard::instance().r_trans;
    TAssembly::Stats stats;
    
    for (const auto &i : r.chromoIDs())
    {
        stats.data[i];
    }
    
    stats.eHist = r.hist();
    stats.iHist = r.hist();
    stats.tHist = r.hist();
    stats.bHist = r.geneHist(ChrT);

    return stats;
}

TAssembly::Stats TAssembly::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_trans;

    assert(!o.ref.empty() && !o.query.empty());
    
    /*
     * 1. Initalize the statistics
     */
    
    TAssembly::Stats stats = init();
    
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

        stats.data[cID].mExonN    = __cmp__.missedExonsN;
        stats.data[cID].mExonR    = __cmp__.missedExonsR;
        stats.data[cID].mExonP    = __cmp__.missedExonsP / 100.0;
        stats.data[cID].mIntronN  = __cmp__.missedIntronsN;
        stats.data[cID].mIntronR  = __cmp__.missedIntronsR;
        stats.data[cID].mIntronP  = __cmp__.missedIntronsP / 100.0;

        stats.data[cID].nExonN    = __cmp__.novelExonsN;
        stats.data[cID].nExonR    = __cmp__.novelExonsR;
        stats.data[cID].nExonP    = __cmp__.novelExonsP / 100.0;
        stats.data[cID].nIntronN  = __cmp__.novelIntronsN;
        stats.data[cID].nIntronR  = __cmp__.novelIntronsR;
        stats.data[cID].nIntronP  = __cmp__.novelIntronsP / 100.0;
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

static void writeChrSummary(const FileName &file, const TAssembly::Stats &stats, const ChromoID &cID, std::shared_ptr<Writer> writer)
{
    const auto &r = Standard::instance().r_trans;
    const auto data = stats.data.at(cID);

    writer->open(file);
    writer->write((boost::format(summary()) % file
                                            % stats.n_expT
                                            % stats.n_chrT
                                            % r.data().size()
                                            % r.countIntrons(cID)
                                            % data.eSN              // 6
                                            % data.eFSN
                                            % data.eSP
                                            % data.eFSP
                                            % (stats.eLimit.id.empty() ? "-" : std::to_string(stats.eLimit.abund))
                                            % stats.eLimit.id
                                            % data.iSN              // 12
                                            % data.iFSN
                                            % data.iSP
                                            % data.iFSP
                                            % (stats.iLimit.id.empty() ? "-" : std::to_string(stats.iLimit.abund))
                                            % stats.iLimit.id
                                            % data.bSN               // 18
                                            % data.bSP
                                            % (stats.bLimit.id.empty() ? "-" : std::to_string(stats.bLimit.abund)) // 20
                                            % stats.bLimit.id        // 21
                                            % data.cSN
                                            % data.cFSN
                                            % data.cSP
                                            % data.cFSP
                                            % data.tSN
                                            % data.tFSN
                                            % data.tSP
                                            % data.tFSP              // 29
                                            % data.mExonN            // 30
                                            % data.mExonR
                                            % data.mExonP
                                            % data.mIntronN
                                            % data.mIntronR          // 34
                                            % data.mIntronP
                                            % data.nExonN
                                            % data.nExonR
                                            % data.nExonP
                                            % data.nIntronN
                                            % data.nIntronR          // 40
                                            % data.nIntronP
                   ).str());
    writer->close();
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

static void writeReplicate(const FileName &file, const TAssembly::Stats &stats, const TAssembly::Options &o)
{
    o.info("Generating statistics for: " + file);

    /*
     * Generating summary statistics for each chromosome
     */
    
    for (const auto &i : stats.data)
    {
        writeChrSummary("TransAssembly_summary (" + i.first + ").stats", stats, i.first, o.writer);
    }
}

void TAssembly::report(const FileName &file, const Options &o)
{
    const auto stats = TAssembly::analyze(file, o);

    /*
     * Generating statistics for the replicate
     */
    
    writeReplicate(file, stats, o);
    
    /*
     * Generating overall statistics
     */
    
    o.info("Generating summary statistics");
    //writeSummary("TransAssembly_summary.stats", file, stats, o);

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
        //writeSummary((boost::format("TransAssembly_%1%_summary.stats") % files[i]).str(), files[i], stats[i], o);
    }
    
    /*
     * Generating sequin statistics for each replicate
     */
    
    for (auto i = 0; i < files.size(); i++)
    {
        writeSequins((boost::format("TransAssembly_%1%_sequin.stats") % files[i]).str(), files[i], stats[i], o);
    }
}