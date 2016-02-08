#include "VARQuin/v_discover.hpp"

using namespace Anaquin;

VDiscover::Stats VDiscover::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VDiscover::Stats stats;

    // Initialize the distribution for each sequin
    stats.hist = r.hist();

    auto &chrT = stats.chrT;

    parseVariant(file, o.caller, [&](const VariantMatch &m)
    {
        if (m.query.chrID == ChrT)
        {
            stats.n_chrT++;
            
            if (m.match && m.ref && m.alt)
            {
                stats.hist.at(m.match->id)++;
                
                if (m.query.pval <= o.sign)
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
                    if (m.query.pval <= o.sign)
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
                                                   % i.query.pval
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

static void writeSummary(const FileName &file, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Number of variants called in the synthetic and experimental chromosomes\n"
                         "   ***\n\n"
                         "   Synthetic:  %2% variants\n\n"
                         "   Experiment: %3% variants\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %4%\n\n"
                         "   Synthetic:  %5% SNPs\n"
                         "   Synthetic:  %6% indels\n"
                         "   Synthetic:  %7% variants\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Experiment)\n"
                         "   ***\n\n"
                         "   File: %8%\n\n"
                         "   Experiment:  %9% SNPs\n"
                         "   Experiment:  %10% indels\n"
                         "   Experiment:  %11% variants\n\n"
                         "   ************************************************************\n"
                         "   ***                                                      ***\n"
                         "   ***        Statistics for the synthetic chromosome       ***\n"
                         "   ***                                                      ***\n"
                         "   ************************************************************\n\n"
                         "   ***\n"
                         "   Detected:    %12% SNPs\n"
                         "   Detected:    %13% indels\n"
                         "   Detected:    %14% variants\n\n"
                         "   Signficiance Level: %15%\n\n"
                         "   Filtered:    %16% SNPs\n"
                         "   Filtered:    %17% indels\n"
                         "   Filtered:    %18% variants\n\n"
                         "   True Positives:   %19% SNPS\n"
                         "   True Positives:   %20% indels\n"
                         "   True Positives:   %21% variants\n\n"
                         "   False Positives:  %22% SNPS\n"
                         "   False Positives:  %23% SNPS\n"
                         "   False Positives:  %24% variants\n\n"
                         "   Performance metrics for SNPs\n\n"
                         "   Sensitivity: %25%\n"
                         "   Specificity: %26%\n\n"
                         "   Performance metrics for indels\n\n"
                         "   Sensitivity: %27%\n"
                         "   Specificity: %28%\n\n"
                         "   Overall performance metrics\n\n"
                         "   Sensitivity: %29%\n"
                         "   Specificity: %30%\n\n";

    o.writer->open("VarDiscover_summary.stats");
    
    o.writer->write((boost::format(summary) % file
                                            % stats.chrT.detectTot()
                                            % stats.endo.size()
                                            % o.rChrT()
                                            % r.countSNPs()
                                            % r.countIndels()
                                            % r.countVars()          // 7
                                            % o.rEndo()
                                            % "NA"
                                            % "NA"
                                            % "NA"
                                            % stats.chrT.detectSNP()
                                            % stats.chrT.detectInd()
                                            % stats.chrT.detectTot()
                                            % o.sign                  // 15
                                            % stats.chrT.filterSNP()
                                            % stats.chrT.filterInd()
                                            % stats.chrT.filterTot()
                                            % stats.chrT.tpSNP()
                                            % stats.chrT.tpInd()
                                            % stats.chrT.tpTot()
                                            % stats.chrT.fpSNP()
                                            % stats.chrT.fpInd()
                                            % stats.chrT.fpTot()
                                            % stats.chrT.m_snp.sn()   // 25 
                                            % stats.chrT.m_snp.sp()
                                            % stats.chrT.m_ind.sn()
                                            % stats.chrT.m_ind.sp()
                                            % stats.chrT.m.sn()
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
    
    writeSummary("VarDiscover_summary.stats", stats, o);
    
    /*
     * Generating labels for the variants
     */
    
    writeClass("VarDiscover_labels.csv", stats.chrT, o);

    /*
     * Generating ROC curve
     */
    
    o.writer->open("VarDiscover_ROC.R");
    o.writer->write(RWriter::createROC_V("VarDiscover_labels.csv"));
    o.writer->close();

    /*
     * Generating LODR curve (only if probability is given, for example, not possible with GATK)
     */

    if (o.caller == Caller::VarScan)
    {
        o.writer->open("VarDiscover_LODR.R");
        o.writer->close();
    }

    /*
     * Generating sequin statistics
     */

    writeSeqins("VarDiscover_quins.csv", stats, o);
}