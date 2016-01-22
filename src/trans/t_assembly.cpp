#include <fstream>
#include "data/compare.hpp"
#include "data/experiment.hpp"
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

static FileName createFilters(const FileName &ref, const FileName &query, const ChromoID &cID)
{
    assert(cID == ChrT || cID == Endo);
    
    std::string line;

    const auto tmp = tmpnam(NULL);
    
    std::ofstream out(tmp);

    auto f = [&](const FileName &file, std::ofstream &out)
    {
        ParserGTF::parse(file, [&](const Feature &f, const std::string &l, const ParserProgress &)
        {
            // Remember the tool pools all the endogenous together
            const auto &id = (f.cID == ChrT ? ChrT : Endo);
            
            if (cID == id)
            {
                out << l << std::endl;
            }
        });
    };
    
    // Generate a filtered query
    f(query, out);
    
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
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %10% (%11%)\n"
           "   Specificity: %12% (%13%)\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %14%\n"
           "   Specificity: %15%\n"
           "   -------------------- Intron Chain level --------------------\n\n"
           "   Sensitivity: %16% (%17%)\n"
           "   Specificity: %18% (%19%)\n\n"
           "   -------------------- Transcript level --------------------\n\n"
           "   Sensitivity: %20% (%21%)\n"
           "   Specificity: %22% (%23%)\n\n"
           "   Missing exons: %24%/%25% (%26%)\n"
           "   Missing introns: %27%/%28% (%29%)\n\n"
           "   Novel exons: %30%/%31% (%32%)\n"
           "   Novel introns: %33%/%34% (%35%)\n\n";
}

static TAssembly::Stats init()
{
    TAssembly::Stats stats;
    
    stats.data[ChrT];
    stats.data[Endo];

    const auto &r = Standard::instance().r_trans;
    
    stats.eHist = r.hist();
    stats.iHist = r.hist();
    stats.tHist = r.hist();
    stats.bHist = r.geneHist(ChrT);

    return stats;
}

TAssembly::Stats TAssembly::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_trans;

    // We'll need the reference annotation for comparison
    assert(!o.ref.empty());
    
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

        stats.data[cID].cSN  = std::min(__cmp__.c_sn  / 100.0, 1.0);
        stats.data[cID].cSP  = std::min(__cmp__.c_sp  / 100.0, 1.0);
        stats.data[cID].cFSN = std::min(__cmp__.c_fsn / 100.0, 1.0);
        stats.data[cID].cFSP = std::min(__cmp__.c_fsp / 100.0, 1.0);
        
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
        /*
         * Unfortunately there's no way to separate statistics between synthetic
         * and endogenous chromosomes. We'll have to create a copy of the reference
         * and query for both.
         */

        const auto created = createFilters(file, o.ref, cID);

        o.logInfo("Filtered query: "     + created + " created");
        o.logInfo("Invoking Cuffcompare: " + created);

        if (cuffcompare_main(o.ref.c_str(), created.c_str()))
        {
            throw std::runtime_error("Failed to analyze the given transcript. Please check the file and try again.");
        }
    };

    o.info("Generating filtered transcript");

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
        
        if (f.cID != Standard::chrT)
        {
            stats.n_endo++;
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
     * 4. Collecting statistics
     */
    
    o.info("Collecting statistics");
    
    stats.eLimit = r.limit(stats.eHist);
    stats.tLimit = r.limit(stats.tHist);
    stats.iLimit = r.limit(stats.iHist);
    stats.bLimit = r.limitGene(stats.bHist);
    
    return stats;
}

static void writeSummary(const FileName &file, const FileName &name, const TAssembly::Stats &stats, const TAssembly::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto data = stats.data.at(ChrT);

    o.info("Generating statistics for: " + name);
    
    // Create the directory if haven't
    o.writer->create(name);

    o.writer->open(name + "/TransAssembly_summary.stats");
    o.writer->write((boost::format(summary()) % file
                                              % stats.n_endo
                                              % stats.n_chrT
                                              % r.data().size()
                                              % r.countIntrons(ChrT)
                                              % data.eSN              // 6
                                              % data.eFSN
                                              % data.eSP
                                              % data.eFSP
                                              % data.iSN              // 10
                                              % data.iFSN
                                              % data.iSP
                                              % data.iFSP
                                              % data.bSN               // 14
                                              % data.bSP
                                              % data.cSN               // 16
                                              % data.cFSN
                                              % data.cSP
                                              % data.cFSP
                                              % data.tSN
                                              % data.tFSN
                                              % data.tSP
                                              % data.tFSP              // 23
                                              % data.mExonN            // 24
                                              % data.mExonR
                                              % data.mExonP
                                              % data.mIntronN
                                              % data.mIntronR          // 28
                                              % data.mIntronP
                                              % data.nExonN
                                              % data.nExonR
                                              % data.nExonP
                                              % data.nIntronN
                                              % data.nIntronR          // 34
                                              % data.nIntronP).str());
    o.writer->close();
}

static void writeSequins(const FileName &file, const FileName &name, const TAssembly::Stats &stats, const TAssembly::Options &o)
{
    o.writer->open(name + "/TransAssembly_sequin.stats");
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

void TAssembly::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TAssembly::analyze(files, o);

    /*
     * Generating summary statistics for each replicate
     */
    
    o.info("Generating summary statistics");

    for (auto i = 0; i < files.size(); i++)
    {
        writeSummary(files[i], o.exp->names().at(i), stats[i], o);
    }
    
    /*
     * Generating sequin statistics for each replicate
     */
    
    for (auto i = 0; i < files.size(); i++)
    {
        writeSequins(files[i], o.exp->names().at(i), stats[i], o);
    }
}