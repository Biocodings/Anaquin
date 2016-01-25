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

static Scripts sSummary()
{
    return "Summary for input: %1%\n\n"
           "   ***\n"
           "   *** Proportion of features mapped to the synthetic and experiment\n"
           "   ***\n\n"
           "   Synthetic:  %3% features\n"
           "   Experiment: %2% features\n\n"
           "   ***\n"
           "   *** Reference annotation (Synthetic)\n"
           "   ***\n\n"
           "   File: %4%\n\n"
           "   Synthetic:  %5% exons\n"
           "   Synthetic:  %6% introns\n\n"
           "   ***\n"
           "   *** Reference annotation (Experiment)\n"
           "   ***\n\n"
           "   File: %7%\n\n"
           "   Experiment:  %8% exons\n"
           "   Experiment:  %9% introns\n\n"
           "   *************************************************\n"
           "   ***                                           ***\n"
           "   ***    Statistics for synthetic chromosome    ***\n"
           "   ***                                           ***\n"
           "   *************************************************\n\n"
           "   ***\n"
           "   *** The following statistics are computed for exact and fuzzy.\n"
           "   *** The fuzzy level is 10 nucleotides.\n"
           "   ***\n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %10% (%11%)\n"
           "   Specificity: %12% (%13%)\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %14% (%15%)\n"
           "   Specificity: %16% (%17%)\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %18%\n"
           "   Specificity: %19%\n\n"
           "   -------------------- Intron Chain level --------------------\n\n"
           "   Sensitivity: %20% (%21%)\n"
           "   Specificity: %22% (%23%)\n\n"
           "   -------------------- Transcript level --------------------\n\n"
           "   Sensitivity: %24% (%25%)\n"
           "   Specificity: %26% (%27%)\n\n"
           "   Missing exons:   %28%/%29% (%30%)\n"
           "   Missing introns: %31%/%32% (%33%)\n\n"
           "   Novel exons:     %34%/%35% (%36%)\n"
           "   Novel introns:   %37%/%38% (%39%)\n\n";
}

static Scripts eSummary()
{
    return "   ***************************************\n"
           "   ***                                 ***\n"
           "   ***    Statistics for experiment    ***\n"
           "   ***                                 ***\n"
           "   ***************************************\n\n"
           "   ***\n"
           "   *** The following statistics are computed for exact and fuzzy.\n"
           "   *** The fuzzy level is 10 nucleotides.\n"
           "   ***\n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %1% (%2%)\n"
           "   Specificity: %3% (%4%)\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %5% (%6%)\n"
           "   Specificity: %7% (%8%)\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %9%\n"
           "   Specificity: %10%\n\n"
           "   -------------------- Intron Chain level --------------------\n\n"
           "   Sensitivity: %11% (%12%)\n"
           "   Specificity: %13% (%14%)\n\n"
           "   -------------------- Transcript level --------------------\n\n"
           "   Sensitivity: %15% (%16%)\n"
           "   Specificity: %17% (%18%)\n\n"
           "   Missing exons:   %19%/%20% (%21%)\n"
           "   Missing introns: %22%/%23% (%24%)\n\n"
           "   Novel exons:   %25%/%26% (%27%)\n"
           "   Novel introns: %28%/%29% (%30%)\n\n";
}

static TAssembly::Stats init(const TAssembly::Options &o)
{
    TAssembly::Stats stats;
    
    stats.data[ChrT];
    stats.refs[ChrT] = o.chrT;
    
    // Remember, endogenous is optional
    if (!o.endo.empty())
    {
        stats.data[Endo];
        stats.refs[Endo] = o.endo;
    }

    return stats;
}

TAssembly::Stats TAssembly::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_trans;

    // We'll need the reference annotation for comparison (endogenous is optional)
    assert(!o.chrT.empty());
    
    /*
     * 1. Initalize the statistics
     */
    
    TAssembly::Stats stats = init(o);
    
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
        // Reference annoation
        const auto &ref = stats.refs.at(cID);

        // Filtered query
        const auto qry = createFilters(file, ref, cID);

        o.logInfo("Reference: " + ref);
        o.logInfo("Query: " + qry);
        
        if (cuffcompare_main(ref.c_str(), qry.c_str()))
        {
            throw std::runtime_error("Failed to analyze " + file + ". Please check the file and try again.");
        }
    };

    o.info("Generating filtered transcript");

    std::for_each(stats.data.begin(), stats.data.end(), [&](const std::pair<ChromoID, TAssembly::Stats::Data> &p)
    {
        compareGTF(p.first);
        copyStats(p.first);
    });
    
    ParserGTF::parse(file, [&](const Feature &f, const std::string &, const ParserProgress &p)
    {
        if (f.cID == ChrT)
        {
            stats.n_chrT++;
        }
        else
        {
            stats.n_endo++;
        }
    });
    
    return stats;
}

static void writeSummary(const FileName &file, const FileName &name, const TAssembly::Stats &stats, const TAssembly::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto data = stats.data.at(ChrT);

    o.info("Generating statistics for: " + name);
    
    #define S(x) (x == 1.0 ? "1.00" : std::to_string(x))
    
    o.writer->create(name);
    o.writer->open(name + "/TransAssembly_summary.stats");
    o.writer->write((boost::format(sSummary()) % file
                                               % stats.n_chrT
                                               % stats.n_endo
                                               % o.rChrT
                                               % r.countExons(ChrT)
                                               % r.countIntrons(ChrT)
                                               % (o.rEndo.empty() ? "-"  : o.rEndo)
                                               % (o.rEndo.empty() ? "NA" : std::to_string(r.countExons("chr1")))
                                               % (o.rEndo.empty() ? "NA" : std::to_string(r.countIntrons("chr1")))
                                               % S(data.eSN)            // 10
                                               % S(data.eFSN)
                                               % S(data.eSP)
                                               % S(data.eFSP)
                                               % S(data.iSN)            // 14
                                               % S(data.iFSN)
                                               % S(data.iSP)
                                               % S(data.iFSP)
                                               % S(data.bSN)            // 18
                                               % S(data.bSP)
                                               % S(data.cSN)            // 20
                                               % S(data.cFSN)
                                               % S(data.cSP)
                                               % S(data.cFSP)
                                               % S(data.tSN)
                                               % S(data.tFSN)
                                               % S(data.tSP)
                                               % S(data.tFSP)           // 27
                                               % data.mExonN            // 28
                                               % data.mExonR
                                               % S(data.mExonP)
                                               % data.mIntronN
                                               % data.mIntronR          // 32
                                               % S(data.mIntronP)
                                               % data.nExonN
                                               % data.nExonR
                                               % S(data.nExonP)
                                               % data.nIntronN
                                               % data.nIntronR          // 38
                                               % S(data.nIntronP)).str());
    if (stats.data.count(Endo))
    {
        const auto &data = stats.data.at(Endo);

        o.writer->write((boost::format(eSummary()) % S(data.eSN)        // 1
                                                   % S(data.eFSN)
                                                   % S(data.eSP)
                                                   % S(data.eFSP)
                                                   % S(data.iSN)        // 5
                                                   % S(data.iFSN)
                                                   % S(data.iSP)
                                                   % S(data.iFSP)
                                                   % S(data.bSN)        // 9
                                                   % S(data.bSP)
                                                   % S(data.cSN)        // 11
                                                   % S(data.cFSN)
                                                   % S(data.cSP)
                                                   % S(data.cFSP)
                                                   % S(data.tSN)
                                                   % S(data.tFSN)
                                                   % S(data.tSP)
                                                   % S(data.tFSP)       // 18
                                                   % data.mExonN        // 19
                                                   % data.mExonR
                                                   % S(data.mExonP)
                                                   % data.mIntronN
                                                   % data.mIntronR      // 23
                                                   % S(data.mIntronP)
                                                   % data.nExonN
                                                   % data.nExonR
                                                   % S(data.nExonP)
                                                   % data.nIntronN
                                                   % data.nIntronR      // 29
                                                   % S(data.nIntronP)).str());
    }

    o.writer->close();
}

//static void writeSequins(const FileName &file, const FileName &name, const TAssembly::Stats &stats, const TAssembly::Options &o)
//{
//    o.writer->open(name + "/TransAssembly_sequin.stats");
//    o.writer->write((boost::format("Summary for file: %1%\n") % file).str());
//    
//    auto format = "%1%\t%2%\t%3%\t%4%";
//    o.writer->write((boost::format(format) % "ID" % "Exon" % "Intron" % "Transcript").str());
//    
//    for (const auto &i : stats.eHist)
//    {
//        o.writer->write((boost::format(format) % i.first
//                                               % stats.eHist.at(i.first)
//                                               % stats.iHist.at(i.first)
//                                               % stats.tHist.at(i.first)).str());
//    }
//
//    o.writer->close();
//}

void TAssembly::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TAssembly::analyze(files, o);

    /*
     * Generating summary statistics for each sample
     */
    
    o.info("Generating summary statistics");

    for (auto i = 0; i < files.size(); i++)
    {
        writeSummary(files[i], o.exp->names().at(i), stats[i], o);
    }
    
    /*
     * Generating sequin statistics for each sample
     */
    
    //for (auto i = 0; i < files.size(); i++)
    //{
    //    writeSequins(files[i], o.exp->names().at(i), stats[i], o);
    //}
}