#include <fstream>
#include "data/path.hpp"
#include "data/compare.hpp"
#include "parsers/parser_gtf.hpp"
#include "RnaQuin/r_assembly.hpp"

#define CHECK_AND_SORT(t) { assert(!t.empty()); std::sort(t.begin(), t.end(), [](const Feature& x, const Feature& y) { return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end); }); }

using namespace Anaquin;

// Defined for cuffcompare
Compare __cmp__;

// Defined in resources.cpp
extern Scripts PlotSensitivity();

// Defined in resources.cpp
extern FileName GTFRef();

// Defined for cuffcompare
extern int cuffcompare_main(const char *ref, const char *query);

template <typename F> inline FileName grepGTF(const FileName &file, F f)
{
    assert(!file.empty());

    Line line;
    
    const auto tmp = tmpFile();
    std::ofstream out(tmp);
    
    auto parse = [&](const FileName &file, std::ofstream &out)
    {
        ParserGTF::parse(file, [&](const Feature &x, const std::string &l, const ParserProgress &)
        {
            if (f(x))
            {
                out << l << std::endl;
            }
        });
    };
    
    parse(file, out);
    
    out.close();
    
    return tmp;
}

// Create a tempatory file for a condition (synthetic or genome)
template <typename F> FileName createGTF(const FileName &file, F f)
{
    Line line;

    const auto tmp = tmpFile();
    std::ofstream out(tmp);

    auto x = [&](const FileName &file, std::ofstream &out)
    {
        ParserGTF::parse(file, [&](const Feature &i, const std::string &l, const ParserProgress &)
        {
            if (f(i.cID))
            {
                out << l << std::endl;
            }
        });
    };
    
    x(file, out);
    out.close();
    
    std::cout << tmp << std::endl;
    return tmp;
}

static FileName createGTFSyn(const FileName &file)
{
    return createGTF(file, [&](const ChrID &cID)
    {
        return Standard::isSynthetic(cID);
    });
}

static FileName createGTFGen(const FileName &file)
{
    return createGTF(file, [&](const ChrID &cID)
    {
        return !Standard::isSynthetic(cID);
    });
}


static RAssembly::Stats init(const RAssembly::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    RAssembly::Stats stats;

    for (const auto &i : r.histGene())
    {
        stats.data[i.first];
    }

    return stats;
}

RAssembly::Stats RAssembly::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_trans;

    auto stats = init(o);

    /*
     * Filtering transcripts
     */
    
    auto copyStats = [&](const ChrID &cID)
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
    
    auto compareGTF = [&](const ChrID &cID, const FileName &ref, const FileName &qry)
    {
        o.logInfo("Reference: " + ref);
        o.logInfo("Query: " + qry);

        #define CUFFCOMPARE(x, y) { if (cuffcompare_main(x.c_str(), y.c_str())) { throw std::runtime_error("Failed to analyze " + file + ". Please check the file and try again."); } }

        // Only required for sensitivity at individual sequins...
        if (Standard::isSynthetic(cID))
        {
            /*
             * Calculating sensitivty for each sequin. Unfortunately, there is no simpler way
             * than Cuffcompare. Cuffcompare doesn't show sensitivity for each isoforms,
             * thus we'll need to generate a GTF file for each sequin.
             */

            const auto hist = r.hist();
            
            for (const auto &i : hist)
            {
                /*
                 * Generate a new GTF solely for the sequin, which will be the reference.
                 */
                
                const auto tmp = grepGTF(ref, [&](const Feature &f)
                {
                    return f.tID == i.first;
                });
                
                o.logInfo("Analyzing: " + i.first);
                
                // Compare only the sequin against the reference
                CUFFCOMPARE(tmp, qry);
                
                stats.tSPs[i.first] = __cmp__.b_sn;
            }
        }

        // Compare everything about the chromosome against the reference
        CUFFCOMPARE(ref, qry);
    };

    o.info("Analyzing transcripts");

    /*
     * Comparing for the synthetic chromosome
     */
    
    compareGTF(ChrT, createGTFSyn(GTFRef()), createGTFSyn(file));
    copyStats(ChrT);
    
    /*
     * Comparing for the genome
     */
    
    if (stats.data.size() > 1)
    {
        compareGTF(Geno, createGTFGen(GTFRef()), createGTFGen(file));
        copyStats(Geno);
    }

    /*
     * Counting exons and transcripts (simply for reporting)
     */
    
    ParserGTF::parse(file, [&](const Feature &f, const std::string &, const ParserProgress &p)
    {
        switch (f.type)
        {
            case Exon:
            {
                if (Standard::isSynthetic(f.cID)) { stats.sExons++; } else { stats.gExons++; } break;
            }

            case Transcript:
            {
                if (Standard::isSynthetic(f.cID)) { stats.sTrans++; } else { stats.gTrans++; } break;
            }

            default: { break; }
        }
    });
    
    return stats;
}

static void generateQuins(const FileName &file, const RAssembly::Stats &stats, const RAssembly::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto format = "%1%\t%2%\t%3%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Seq"
                                           % "Expected"
                                           % "Measured").str());

    for (const auto &i : stats.tSPs)
    {
        o.writer->write((boost::format(format) % i.first
                                               % r.match(i.first)->concent()
                                               % (stats.tSPs.at(i.first) / 100.0)).str());
    }
    
    o.writer->close();
}

static void generateSummary(const FileName &file, const RAssembly::Stats &stats, const RAssembly::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    const auto hasGeno = stats.data.size() > 1;
    
    const auto sData = stats.data.at(ChrT);
    const auto gData = hasGeno ? stats.data.at(Geno) : RAssembly::Stats::Data();

    #define S(x) (x == 1.0 ? "1.00" : std::to_string(x))
    
    const auto format = "-------RnaAssembly Summary Statistics\n\n"
                        "       User assembly file: %1%\n"
                        "       Reference annotation file: %2%\n\n"
                        "-------User Gene Assemblies (Synthetic)\n\n"
                        "       Synthetic: %3% exons\n"
                        "       Synthetic: %4% introns\n"
                        "       Synthetic: %5% isoforms\n"
                        "       Synthetic: %6% genes\n\n"
                        "-------User Gene Assemblies (Genome)\n\n"
                        "       Genome: %7% exons\n"
                        "       Genome: %8% introns\n"
                        "       Genome: %9% isoforms\n"
                        "       Genome: %10% genes\n\n"
                        "-------Reference Gene Annotations (Synthetic)\n\n"
                        "       Synthetic: %11% exons\n"
                        "       Synthetic: %12% introns\n"
                        "       Synthetic: %13% isoforms\n"
                        "       Synthetic: %14% genes\n\n"
                        "-------Reference Gene Annotations (Genome)\n\n"
                        "       Genome: %15% exons\n"
                        "       Genome: %16% introns\n"
                        "       Genome: %17% isoforms\n"
                        "       Genome: %18% genes\n\n"
                        "-------Comparison of assembly to annotations (Synthetic)\n\n"
                        "       *Exon level\n"
                        "       Sensitivity: %19%\n"
                        "       Specificity: %20%\n\n"
                        "       *Intron\n"
                        "       Sensitivity: %21%\n"
                        "       Specificity: %22%\n\n"
                        "       *Base level\n"
                        "       Sensitivity: %23%\n"
                        "       Specificity: %24%\n\n"
                        "       *Intron Chain\n"
                        "       Sensitivity: %25%\n"
                        "       Specificity: %26%\n\n"
                        "       *Transcript level\n"
                        "       Sensitivity: %27%\n"
                        "       Specificity: %28%\n\n"
                        "       Missing exons: %29%\n"
                        "       Missing introns: %30%\n\n"
                        "       Novel exons:  %31%\n"
                        "       Novel introns: %32%\n\n"
                        "-------Comparison of assembly to annotations (Genome)\n\n"
                        "       *Exon level\n"
                        "       Sensitivity: %33%\n"
                        "       Specificity: %34%\n\n"
                        "       *Intron\n"
                        "       Sensitivity: %35%\n"
                        "       Specificity: %36%\n\n"
                        "       *Base level\n"
                        "       Sensitivity: %37%\n"
                        "       Specificity: %38%\n\n"
                        "       *Intron Chain\n"
                        "       Sensitivity: %39%\n"
                        "       Specificity: %40%\n\n"
                        "       *Transcript level\n"
                        "       Sensitivity: %41%\n"
                        "       Specificity: %42%\n\n"
                        "       Missing exons: %43%\n"
                        "       Missing introns: %44%\n\n"
                        "       Novel exons:  %45%\n"
                        "       Novel introns: %46%";
    
    o.generate("RnaAssembly_summary.stats");
    o.writer->open("RnaAssembly_summary.stats");
    o.writer->write((boost::format(format) % file
                                           % GTFRef()
                                           % stats.sExons
                                           % stats.sIntrons
                                           % stats.sTrans
                                           % stats.sGenes
                                           % stats.gExons
                                           % stats.gIntrons
                                           % stats.gTrans
                                           % stats.gGenes      // 10
                                           % r.countExonSyn()  // 11
                                           % r.countIntrSyn()  // 12
                                           % r.countTransSyn() // 13
                                           % r.countGeneSyn()  // 14
                                           % r.countExonSyn()  // 15
                                           % r.countIntrSyn()  // 16
                                           % r.countTransSyn() // 17
                                           % r.countGeneSyn()  // 18
                                           % S(sData.eSN)      // 19
                                           % S(sData.eSP)      // 20
                                           % S(sData.iSN)      // 21
                                           % S(sData.iSP)      // 22
                                           % S(sData.bSN)      // 23
                                           % S(sData.bSP)      // 24
                                           % S(sData.cSN)      // 25
                                           % S(sData.cSP)      // 26
                                           % S(sData.tSN)      // 27
                                           % S(sData.tSP)      // 28
                                           % sData.mExonN      // 29
                                           % sData.mIntronN    // 30
                                           % sData.nExonN      // 31
                                           % sData.nIntronN    // 32
                                           % (hasGeno ? S(gData.eSN)      : "-") // 33
                                           % (hasGeno ? S(gData.eSP)      : "-") // 34
                                           % (hasGeno ? S(gData.iSN)      : "-") // 35
                                           % (hasGeno ? S(gData.iSP)      : "_") // 36
                                           % (hasGeno ? S(gData.bSN)      : "-") // 37
                                           % (hasGeno ? S(gData.bSP)      : "-") // 38
                                           % (hasGeno ? S(gData.cSN)      : "-") // 39
                                           % (hasGeno ? S(gData.cSP)      : "-") // 40
                                           % (hasGeno ? S(gData.tSN)      : "-") // 41
                                           % (hasGeno ? S(gData.tSN)      : "-") // 42
                                           % (hasGeno ? S(gData.mExonN)   : "-") // 43
                                           % (hasGeno ? S(gData.mIntronN) : "-") // 44
                                           % (hasGeno ? S(gData.nExonN)   : "-") // 45
                                           % (hasGeno ? S(gData.nIntronN) : "-")).str());
    o.writer->close();
}

void RAssembly::report(const FileName &file, const Options &o)
{
    const auto stats = RAssembly::analyze(file, o);

    /*
     * Generating RnaAssembly_summary.stats
     */
    
    generateSummary(file, stats, o);
    
    /*
     * Generating RnaAssembly_quins.stats
     */

    generateQuins("RnaAssembly_quins.csv", stats, o);
    
    /*
     * Generating RnaAssembly_assembly.R
     */
    
    o.generate("RnaAssembly_assembly.R");
    o.writer->open("RnaAssembly_assembly.R");
    o.writer->write(RWriter::createSensitivity("RnaAssembly_quins.csv",
                                               "Assembly Detection",
                                               "Input Concentration (log2)",
                                               "Sensitivity",
                                               "Expected",
                                               "Measured",
                                               true));
    o.writer->close();
    
    /*
     * Generating RnaAssembly_report.pdf
     */
    
    o.report->open("RnaAssembly_report.pdf");
    o.report->addTitle("RnaAssembly");
    o.report->addFile("RnaAssembly_summary.stats");
    o.report->addFile("RnaAssembly_quins.csv");
    o.report->addFile("RnaAssembly_assembly.R");
}