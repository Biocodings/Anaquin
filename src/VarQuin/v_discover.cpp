#include "VarQuin/v_discover.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotLODR_V();

// Defined in resources.cpp
extern Scripts PlotROC_V();

VDiscover::Stats VDiscover::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VDiscover::Stats stats;

    // Initialize the distribution for each sequin
    stats.hist = r.hist();

    auto &chrT = stats.chrT;

    parseVariant(file, o.caller, [&](const VariantMatch &m)
    {
        if (m.query.cID == ChrT)
        {
            stats.n_chrT++;
            
            if (m.match && m.ref && m.alt)
            {
                stats.hist.at(m.match->id)++;
                
                if (m.query.p <= o.sign)
                {
                    chrT.tps.push_back(m);
                }
                else
                {
                    chrT.tns.push_back(m);
                }
            }
            else
            {
                if (!m.seq)
                {
                    o.warn("Variant [" + std::to_string(m.query.l.start) + "] not aligned with any of the sequins. Ignored.");
                }
                else
                {
                    if (m.query.p <= o.sign)
                    {
                        chrT.fps.push_back(m);
                    }
                    else
                    {
                        chrT.fns.push_back(m);
                    }
                }
            }
        }
        else
        {
            stats.n_endo++;
            stats.endo.push_back(m.query);
        }
    });

    for (const auto &i : chrT.tps)
    {
        chrT.m.tp()++;
        
        switch (i.query.type())
        {
            case Mutation::SNP:       { chrT.m_snp.tp()++; break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { chrT.m_ind.tp()++; break; }
        }
    }

    for (const auto &i : chrT.fps)
    {
        chrT.m.fp()++;
        
        switch (i.query.type())
        {
            case Mutation::SNP:       { chrT.m_snp.fp()++; break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { chrT.m_ind.fp()++; break; }
        }
    }

    for (const auto &i : chrT.tns)
    {
        chrT.m.tn()++;

        switch (i.query.type())
        {
            case Mutation::SNP:       { chrT.m_snp.tn()++; break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { chrT.m_ind.tn()++; break; }
        }
    }

    for (const auto &i : chrT.fns)
    {
        chrT.m.fn()++;
        
        switch (i.query.type())
        {
            case Mutation::SNP:       { chrT.m_snp.fn()++; break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { chrT.m_ind.fn()++; break; }
        }
    }
    
    chrT.m_snp.nq() = chrT.detectSNP();
    chrT.m_snp.nr() = r.countSNPs();

    chrT.m_ind.nq() = chrT.detectInd();
    chrT.m_ind.nr() = r.countIndels();

    chrT.m.nq() = stats.chrT.m_snp.nq() + stats.chrT.m_ind.nq();
    chrT.m.nr() = stats.chrT.m_snp.nr() + stats.chrT.m_ind.nr();
    
    return stats;
}

static void writeClass(const FileName &file,
                       const VDiscover::Stats::ChrTStats &stats,
                       const VDiscover::Options &o)
{
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Sequin"
                                           % "Position"
                                           % "Label"
                                           % "PValue"
                                           % "RefRead"
                                           % "VarRead"
                                           % "ERatio"
                                           % "EAlleleF"
                                           % "Type").str());

    auto f = [&](const std::vector<VDiscover::Stats::ChrTData> &x, const std::string &label)
    {
        for (const auto &i : x)
        {
            std::string type;
            
            switch (i.query.type())
            {
                case Mutation::SNP:       { type = "SNP";   break; }
                case Mutation::Deletion:
                case Mutation::Insertion: { type = "Indel"; break; }
            }
            
            o.writer->write((boost::format(format) % i.seq->id
                                                   % i.query.l.start
                                                   % label
                                                   % i.query.p
                                                   % i.query.readR
                                                   % i.query.readV
                                                   % i.eFold
                                                   % i.eAllFreq
                                                   % type).str());
        }
    };

    f(stats.tps, "TP");
    f(stats.fps, "FP");
    f(stats.tns, "TN");
    f(stats.fns, "FN");

    o.writer->close();
}

static void writeSeqins(const FileName &file, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    o.writer->open(file);
    o.writer->write((boost::format("Summary for input: %1%\n") % file).str());
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";
    o.writer->write((boost::format(format) % "Sequin"
                                           % "RefSNP"
                                           % "RefIndel"
                                           % "RefTotal"
                                           % "FoundSNP"
                                           % "FoundIndel"
                                           % "FoundTotal"
                                           % "SN (SNP)"
                                           % "SP (SNP"
                                           % "SN (Indel)"
                                           % "SP (Indel)"
                                           % "SN"
                                           % "SP").str());
    o.writer->close();
}

static void writeSummary(const FileName &file, const FileName &src, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Number of variants called in the synthetic and experimental chromosomes\n"
                         "   ***\n\n"
                         "   Synthetic:  %2% variants\n"
                         "   Experiment: %3% variants\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %4%\n\n"
                         "   Synthetic:  %5% SNPs\n"
                         "   Synthetic:  %6% indels\n"
                         "   Synthetic:  %7% variants\n\n"
                         "   ************************************************************\n"
                         "   ***                                                      ***\n"
                         "   ***        Statistics for the synthetic chromosome       ***\n"
                         "   ***                                                      ***\n"
                         "   ************************************************************\n\n"
                         "   Detected:    %8% SNPs\n"
                         "   Detected:    %9% indels\n"
                         "   Detected:    %10% variants\n\n"
                         "   Signficiance Level: %11%\n\n"
                         "   Filtered:    %12% SNPs\n"
                         "   Filtered:    %13% indels\n"
                         "   Filtered:    %14% variants\n\n"
                         "   True Positives:   %15% SNPS\n"
                         "   True Positives:   %16% indels\n"
                         "   True Positives:   %17% variants\n\n"
                         "   False Positives:  %18% SNPS\n"
                         "   False Positives:  %19% SNPS\n"
                         "   False Positives:  %20% variants\n\n"
                         "   ***\n"
                         "   *** Performance metrics (Overall)\n"
                         "   ***\n\n"
                         "   Sensitivity: %21%\n"
                         "   Specificity: %22%\n\n"
                         "   ***\n"
                         "   *** Performance metrics (SNP)\n"
                         "   ***\n\n"
                         "   Sensitivity: %23%\n"
                         "   Specificity: %24%\n\n"
                         "   ***\n"
                         "   *** Performance metrics (Indel)\n"
                         "   ***\n\n"
                         "   Sensitivity: %25%\n"
                         "   Specificity: %26%\n\n";

    o.writer->open("VarDiscover_summary.stats");    
    o.writer->write((boost::format(summary) % src
                                            % stats.chrT.detectTot()
                                            % stats.endo.size()
                                            % o.rChrT
                                            % r.countSNPs()
                                            % r.countIndels()
                                            % r.countVars()           // 7
                                            % stats.chrT.detectSNP()
                                            % stats.chrT.detectInd()
                                            % stats.chrT.detectTot()
                                            % o.sign                  // 11
                                            % stats.chrT.filterSNP()
                                            % stats.chrT.filterInd()
                                            % stats.chrT.filterTot()
                                            % stats.chrT.tpSNP()
                                            % stats.chrT.tpInd()
                                            % stats.chrT.tpTot()
                                            % stats.chrT.fpSNP()
                                            % stats.chrT.fpInd()
                                            % stats.chrT.fpTot()
                                            % stats.chrT.m_snp.sn()   // 21
                                            % stats.chrT.m_snp.sp()   // 22
                                            % stats.chrT.m_ind.sn()   // 23
                                            % stats.chrT.m_ind.sp()   // 24
                                            % stats.chrT.m.sn()       // 25
                                            % stats.chrT.m.sp()).str());
    o.writer->close();
}

void VDiscover::report(const FileName &file, const Options &o)
{
    o.logInfo("Significance level: " + std::to_string(o.sign));

    const auto stats = analyze(file, o);

    o.logInfo("Number of true positives:  " + std::to_string(stats.chrT.tps.size()));
    o.logInfo("Number of false positives: " + std::to_string(stats.chrT.fps.size()));
    o.logInfo("Number of true negatives: "  + std::to_string(stats.chrT.tns.size()));
    o.logInfo("Number of false negaitves: " + std::to_string(stats.chrT.fns.size()));

    o.info("Generating statistics");

    /*
     * Generate summary statistics
     */
    
    writeSummary("VarDiscover_summary.stats", file, stats, o);

    /*
     * Generating classified statistics for the variants
     */
    
    writeClass("VarDiscover_labels.csv", stats.chrT, o);

    /*
     * Generating ROC curve
     */
    
    o.writer->open("VarDiscover_ROC.R");
    o.writer->write(RWriter::createScript("VarDiscover_labels.csv", PlotROC_V()));
    o.writer->close();

    /*
     * Generating LODR curve (only if probability is given, for instance, not possible with GATK)
     */

    if (o.caller == Caller::VarScan)
    {
        o.writer->open("VarDiscover_LODR.R");
        o.writer->write(RWriter::createScript("VarDiscover_labels.csv", PlotLODR_V()));
        o.writer->close();
    }

    /*
     * Generating sequin statistics
     */

    writeSeqins("VarDiscover_quins.csv", stats, o);
}