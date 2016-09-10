#include <thread>
#include <fstream>
#include "data/file.hpp"
#include "data/compare.hpp"
#include "tools/gtf_data.hpp"
#include "RnaQuin/r_assembly.hpp"

using namespace Anaquin;

// Defined for cuffcompare
Compare __cmp__;

static FileName __QForSyn__;
static FileName __QForGen__;
static FileName __RForSyn__;
static FileName __RForGen__;

static GTFData __RData__;
static RAssembly::Stats *__Stats__;

static bool __hasGen__;

// Defined in resources.cpp
extern FileName GTFRef();

// Defined for Cuffcompare
extern int cuffcompare_main(const char *ref, const char *query);

static std::string exec(const char* cmd)
{
    char buffer[128];
    std::string result = "";
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    
    return result;
}

static FileName grepGTF(const FileName &file, const std::string &key, bool shouldGeneTrans = true)
{
    const auto tmp = tmpFile();

    std::string cmd;
    
    if (shouldGeneTrans)
    {
        cmd = "grep " + key + " " + file + " > " + tmp;
    }
    else
    {
        cmd = "grep " + key + " " + file + " | grep -v gene_type | grep -v transcript_type > " + tmp;
    }
    
    const auto msg = exec(cmd.c_str());
    
    if (!msg.empty())
    {
        A_THROW("Failed to execute: " + msg);
    }

    return tmp;
}

static FileName grepVGTF(const FileName &file, const std::string &key)
{
    const auto tmp = tmpFile();
    const auto cmd = "grep -v " + key + " " + file + " > " + tmp;
    const auto msg = exec(cmd.c_str());
    
    if (!msg.empty())
    {
        A_THROW("Failed to execute: " + msg);
    }
    
    return tmp;
}

static FileName createQGTFSyn(const FileName &file)
{
    return (__QForSyn__ = grepGTF(file, ChrIS));
}

static FileName createQGTFGen(const FileName &file)
{
    return (__QForGen__ = grepVGTF(file, ChrIS));
}

static FileName createRGTFSyn(const FileName &file)
{
    return (__RForSyn__ = grepGTF(file, ChrIS, false));
}

static FileName createRGTFGen(const FileName &file)
{
    return (__RForGen__ = grepVGTF(file, ChrIS));
}

static void readQueryGTF(const FileName &file)
{
    const auto gs = gtfData(Reader(file));
    
    __Stats__->sExons = gs.countUExonSyn();
    __Stats__->sIntrs = gs.countUIntrSyn();
    __Stats__->sTrans = gs.countTransSyn();
    __Stats__->sGenes = gs.countGeneSyn();
    __Stats__->gExons = gs.countUExonGen();
    __Stats__->gIntrs = gs.countUIntrGen();
    __Stats__->gTrans = gs.countTransGen();
    __Stats__->gGenes = gs.countGeneGen();
}

static void readRefGTF(const FileName &file)
{
    __RData__ = gtfData(Reader(file));
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
    __Stats__ = &stats;

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
                // Generate a new reference GTF solely for the sequins
                const auto tmp = grepGTF(ref, i.first);
                
                o.analyze(i.first);
                
                // Compare only the sequin against the reference
                CUFFCOMPARE(tmp, qry);
                
                stats.tSPs[i.first] = __cmp__.b_sn;
            }
        }

        // Compare everything about the chromosome against the reference
        CUFFCOMPARE(ref, qry);

        o.logInfo("Compare complated");
    };

    std::thread t1(readQueryGTF, file);
    std::thread t2(readRefGTF, GTFRef());

    /*
     * Create a new synthetic-only GTF file, which we'll use to estimate sequin sensitivity
     */

    o.info("Analyzing transcripts");
    o.info("Creating temporary transcripts");
    
    std::thread t3(createQGTFSyn, file);
    std::thread t4(createQGTFGen, file);
    std::thread t5(createRGTFSyn, GTFRef());
    std::thread t6(createRGTFGen, GTFRef());

    t3.join();
    t4.join();
    t5.join();
    t6.join();
    
    o.logInfo(__QForSyn__);
    o.logInfo(__QForGen__);
    o.logInfo(__RForSyn__);
    o.logInfo(__RForGen__);
    
    A_CHECK(!__QForSyn__.empty(), "Error in __QForSyn__");
    A_CHECK(!__QForGen__.empty(), "Error in __QForGen__");
    A_CHECK(!__RForSyn__.empty(), "Error in __RForSyn__");
    A_CHECK(!__RForGen__.empty(), "Error in __RForGen__");
    
    o.info("Temporary transcripts generated");

    /*
     * Comparing for the synthetic
     */

    o.info("Generating for the synthetic");
    compareGTF(ChrIS, __RForSyn__, __QForSyn__);
    copyStats(ChrIS);
    
    /*
     * Comparing for the genome
     */
    
    if ((__hasGen__ = !isEmpty(__QForGen__)))
    {
        o.analyze("Genome");
        compareGTF(Geno, __RForGen__, __QForGen__);
        copyStats(Geno);
    }
    
    o.info("Waiting for worker threads to complete");

    t1.join();
    t2.join();
    
    A_CHECK(stats.sExons, "stats.sExons");

    return stats;
}

static void generateQuins(const FileName &file, const RAssembly::Stats &stats, const RAssembly::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto format = "%1%\t%2%\t%3%\t%4%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Length"
                                           % "InputConcent"
                                           % "Sn").str());

    for (const auto &i : stats.tSPs)
    {
        const auto m = r.findTrans(ChrIS, i.first);
        o.writer->write((boost::format(format) % i.first
                                               % m->l.length()
                                               % r.match(i.first)->concent()
                                               % (stats.tSPs.at(i.first) / 100.0)).str());
    }
    
    o.writer->close();
}   

static void generateSummary(const FileName &file, const RAssembly::Stats &stats, const RAssembly::Options &o)
{
    const auto &r = __RData__;

    const auto hasGen = __hasGen__;
    const auto sData  = stats.data.at(ChrIS);
    const auto gData  = hasGen ? stats.data.at(Geno) : RAssembly::Stats::Data();

    #define C(x) (std::to_string(x))
    #define S(x) (x == 1.0 ? "1.00" : std::to_string(x))
    
    const auto format = "-------RnaAssembly Summary Statistics\n\n"
                        "       User assembly file: %1%\n"
                        "       Reference annotation file: %2%\n\n"
                        "-------Reference Gene Annotations (Synthetic)\n\n"
                        "       Synthetic: %3% exons\n"
                        "       Synthetic: %4% introns\n"
                        "       Synthetic: %5% isoforms\n"
                        "       Synthetic: %6% genes\n\n"
                        "-------Reference Gene Annotations (Genome)\n\n"
                        "       Genome: %7% exons\n"
                        "       Genome: %8% introns\n"
                        "       Genome: %9% isoforms\n"
                        "       Genome: %10% genes\n\n"
                        "-------User Assemblies (Synthetic)\n\n"
                        "       Synthetic: %11% exons\n"
                        "       Synthetic: %12% introns\n"
                        "       Synthetic: %13% isoforms\n"
                        "       Synthetic: %14% genes\n\n"
                        "-------User Assemblies (Genome)\n\n"
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
                        "       Missing exons:   %29%\n"
                        "       Missing introns: %30%\n\n"
                        "       Novel exons:   %31%\n"
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
                        "       Missing exons:   %43%\n"
                        "       Missing introns: %44%\n\n"
                        "       Novel exons:   %45%\n"
                        "       Novel introns: %46%";
    
    o.generate("RnaAssembly_summary.stats");
    o.writer->open("RnaAssembly_summary.stats");
    o.writer->write((boost::format(format) % file              // 1
                                           % GTFRef()          // 2
                                           % r.countUExonSyn() // 3
                                           % r.countUIntrSyn() // 4
                                           % r.countTransSyn() // 5
                                           % r.countGeneSyn()  // 6
                                           % r.countUExonGen() // 7
                                           % r.countUIntrGen() // 8
                                           % r.countTransGen() // 9
                                           % r.countGeneGen()  // 10
                                           % stats.sExons      // 11
                                           % stats.sIntrs      // 12
                                           % stats.sTrans      // 13
                                           % stats.sGenes      // 14
                                           % stats.gExons      // 15
                                           % stats.gIntrs      // 16
                                           % stats.gTrans      // 17
                                           % stats.gGenes      // 18
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
                                           % (hasGen ? S(gData.eSN)      : "-") // 33
                                           % (hasGen ? S(gData.eSP)      : "-") // 34
                                           % (hasGen ? S(gData.iSN)      : "-") // 35
                                           % (hasGen ? S(gData.iSP)      : "_") // 36
                                           % (hasGen ? S(gData.bSN)      : "-") // 37
                                           % (hasGen ? S(gData.bSP)      : "-") // 38
                                           % (hasGen ? S(gData.cSN)      : "-") // 39
                                           % (hasGen ? S(gData.cSP)      : "-") // 40
                                           % (hasGen ? S(gData.tSN)      : "-") // 41
                                           % (hasGen ? S(gData.tSN)      : "-") // 42
                                           % (hasGen ? C(gData.mExonN)   : "-") // 43
                                           % (hasGen ? C(gData.mIntronN) : "-") // 44
                                           % (hasGen ? C(gData.nExonN)   : "-") // 45
                                           % (hasGen ? C(gData.nIntronN) : "-")).str());
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

    generateQuins("RnaAssembly_sequins.csv", stats, o);
    
    /*
     * Generating RnaAssembly_assembly.R
     */
    
    o.generate("RnaAssembly_assembly.R");
    o.writer->open("RnaAssembly_assembly.R");
    o.writer->write(RWriter::createSensitivity("RnaAssembly_sequins.csv",
                                               "Assembly Detection",
                                               "Input Concentration (log2)",
                                               "Sensitivity",
                                               "InputConcent",
                                               "Sn",
                                               true));
    o.writer->close();
    
    /*
     * Generating RnaAssembly_report.pdf
     */
    
    o.report->open("RnaAssembly_report.pdf");
    o.report->addTitle("RnaAssembly");
    o.report->addFile("RnaAssembly_summary.stats");
    o.report->addFile("RnaAssembly_sequins.csv");
    o.report->addFile("RnaAssembly_assembly.R");
}