#include "VarQuin/v_discover.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVProb();

// Defined in resources.cpp
extern Scripts PlotVROC();

VDiscover::Stats VDiscover::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VDiscover::Stats stats;
    stats.hist = r.varHist();

    parseVariant(file, o.soft, [&](const VariantMatch &m)
    {
        if (m.query.cID == ChrT)
        {
            stats.n_chrT++;
            
            const auto p = isnan(m.query.p) ? 0.0 : m.query.p;
            
            if (m.match && m.ref && m.alt)
            {
                const auto key = var2hash(m.match->id, m.match->type(), m.match->l);
                stats.hist.at(key)++;
                
                if (p <= o.sign)
                {
                    stats.chrT.tps.push_back(m);
                }
                else
                {
                    stats.chrT.fns.push_back(m);
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
                    if (p <= o.sign)
                    {
                        stats.chrT.fps.push_back(m);
                    }
                    else
                    {
                        stats.chrT.tns.push_back(m);
                    }
                }
            }
        }
        else
        {
            stats.n_geno++;
            stats.geno.push_back(m.query);
        }
    });

    for (const auto &i : stats.chrT.tps)
    {
        stats.chrT.m.tp()++;
        
        switch (i.query.type())
        {
            case Mutation::SNP:       { stats.chrT.m_snp.tp()++; break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { stats.chrT.m_ind.tp()++; break; }
        }
    }

    for (const auto &i : stats.chrT.fps)
    {
        stats.chrT.m.fp()++;
        
        switch (i.query.type())
        {
            case Mutation::SNP:       { stats.chrT.m_snp.fp()++; break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { stats.chrT.m_ind.fp()++; break; }
        }
    }

    for (const auto &i : stats.chrT.tns)
    {
        stats.chrT.m.tn()++;

        switch (i.query.type())
        {
            case Mutation::SNP:       { stats.chrT.m_snp.tn()++; break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { stats.chrT.m_ind.tn()++; break; }
        }
    }

    for (const auto &i : stats.chrT.fns)
    {
        stats.chrT.m.fn()++;
        
        switch (i.query.type())
        {
            case Mutation::SNP:       { stats.chrT.m_snp.fn()++; break; }
            case Mutation::Deletion:
            case Mutation::Insertion: { stats.chrT.m_ind.fn()++; break; }
        }
    }
    
    stats.chrT.m_snp.nq() = stats.chrT.dSNP();
    stats.chrT.m_snp.nr() = r.countSNPs();

    stats.chrT.m_ind.nq() = stats.chrT.dInd();
    stats.chrT.m_ind.nr() = r.countIndels();

    stats.chrT.m.nq() = stats.chrT.m_snp.nq() + stats.chrT.m_ind.nq();
    stats.chrT.m.nr() = stats.chrT.m_snp.nr() + stats.chrT.m_ind.nr();
    
    return stats;
}

static std::string type2str(Mutation type)
{
    switch (type)
    {
        case Mutation::SNP:       { return "SNP"; }
        case Mutation::Deletion:
        case Mutation::Insertion: { return "Indel"; }
    }
}

static void writeQuins(const FileName &file,
                       const VDiscover::Stats &stats,
                       const VDiscover::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";

    o.writer->open("VarDiscover_quins.stats");
    o.writer->write((boost::format(format) % "seq"
                                           % "pos"
                                           % "label"
                                           % "ref"
                                           % "var"
                                           % "eFold"
                                           % "eAllele"
                                           % "pval"
                                           % "type").str());

    for (const auto &i : stats.hist)
    {
        const auto &r = Standard::instance().r_var;
        
        if (i.second)
        {
            auto f = [&](const std::vector<VDiscover::Stats::ChrTData> &x, const std::string &label)
            {
                for (const auto &j : x)
                {
                    if (i.first == j.seq->key())
                    {
                        o.writer->write((boost::format(format) % j.seq->id
                                                               % j.query.l.start
                                                               % label
                                                               % j.query.readR
                                                               % j.query.readV
                                                               % j.eFold
                                                               % j.eAllFreq
                                                               % (isnan(j.query.p) ? "-" : p2str(j.query.p))
                                                               % type2str(j.query.type())).str());
                        return true;
                    }
                }
                
                return false;
            };
            
            if (!f(stats.chrT.tps, "TP") &&
                !f(stats.chrT.fps, "FP") &&
                !f(stats.chrT.tns, "TN") &&
                !f(stats.chrT.fns, "FN"))
            {
                throw std::runtime_error("Failed to find hash key in writeQuins()");
            }
        }
        else
        {
            const auto m = r.hashVar(i.first);
            assert(m);
            
            o.writer->write((boost::format(format) % m->id
                                                   % m->l.start
                                                   % "FN"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % "NA"
                                                   % type2str(m->type())).str());
        }
    }
    
    o.writer->close();
}

static void writeQuery(const FileName &file, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "seq"
                                           % "pos"
                                           % "label"
                                           % "ref"
                                           % "var"
                                           % "eFold"
                                           % "eAllele"
                                           % "pval"
                                           % "type").str());

    auto f = [&](const std::vector<VDiscover::Stats::ChrTData> &x, const std::string &label)
    {
        for (const auto &i : x)
        {
            o.writer->write((boost::format(format) % i.seq->id
                                                   % i.query.l.start
                                                   % label
                                                   % i.query.readR
                                                   % i.query.readV
                                                   % i.eFold
                                                   % i.eAllFreq
                                                   % (isnan(i.query.p) ? "-" : p2str(i.query.p))
                                                   % type2str(i.query.type())).str());
        }
    };

    f(stats.chrT.tps, "TP");
    f(stats.chrT.fps, "FP");
    f(stats.chrT.tns, "TN");
    f(stats.chrT.fns, "FN");

    o.writer->close();
}

static void writeSummary(const FileName &file, const FileName &src, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    // Can specificity (requires p-value) shown?
    const auto noSP = o.soft != VDiscover::Software::VarScan;
    
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Number of variants called in the synthetic and genome\n"
                         "   ***\n\n"
                         "   Synthetic: %2% variants\n"
                         "   Genome:    %3% variants\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %4%\n\n"
                         "   Synthetic:  %5% SNPs\n"
                         "   Synthetic:  %6% indels\n"
                         "   Synthetic:  %7% variants\n\n"
                         "   ***                                         \n"
                         "   *** Statistics for the synthetic chromosome \n"
                         "   ***                                         \n\n"
                         "   Detected:    %8% SNPs\n"
                         "   Detected:    %9% indels\n"
                         "   Detected:    %10% variants\n\n"
                         "   Significance: %11%\n\n"
                         "   Significant: %12% SNPs\n"
                         "   Significant: %13% indels\n"
                         "   Significant: %14% variants\n\n"
                         "   True Positive:  %15% SNPS\n"
                         "   True Positive:  %16% indels\n"
                         "   True Positive:  %17% variants\n\n"
                         "   False Positive: %18% SNPS\n"
                         "   False Positive: %19% SNPS\n"
                         "   False Positive: %20% variants\n\n"
                         "   ***\n"
                         "   *** Performance metrics (Overall)\n"
                         "   ***\n\n"
                         "   Sensitivity: %21$.2f\n"
                         "   Specificity: %22$.2f\n"
                         "   Precision:   %23$.2f\n\n"
                         "   ***\n"
                         "   *** Performance metrics (SNP)\n"
                         "   ***\n\n"
                         "   Sensitivity: %24$.2f\n"
                         "   Specificity: %25$.2f\n"
                         "   Precision:   %26$.2f\n\n"
                         "   ***\n"
                         "   *** Performance metrics (Indel)\n"
                         "   ***\n\n"
                         "   Sensitivity: %27$.2f\n"
                         "   Specificity: %28$.2f\n"
                         "   Precision:   %29$.2f\n\n";

    o.info("Generating " + file);
    o.writer->open("VarDiscover_summary.stats");
    o.writer->write((boost::format(summary) % src
                                            % stats.chrT.dTot()
                                            % stats.geno.size()
                                            % o.rChrT
                                            % r.countSNPs()
                                            % r.countIndels()
                                            % r.countVars()           // 7
                                            % stats.chrT.dSNP()
                                            % stats.chrT.dInd()
                                            % stats.chrT.dTot()
                                            % o.sign                  // 11
                                            % stats.chrT.sSNP()
                                            % stats.chrT.sInd()
                                            % stats.chrT.sTot()
                                            % stats.chrT.tpSNP()
                                            % stats.chrT.tpInd()
                                            % stats.chrT.tpTot()
                                            % stats.chrT.fpSNP()
                                            % stats.chrT.fpInd()
                                            % stats.chrT.fpTot()
                                            % stats.chrT.m.sn()                              // 21
                                            % (noSP ? "-" : toString(stats.chrT.m.sp()))     // 22
                                            % stats.chrT.m.pc()                              // 23
                                            % stats.chrT.m_snp.sn()                          // 24
                                            % (noSP ? "-" : toString(stats.chrT.m_snp.sp())) // 25
                                            % stats.chrT.m_snp.pc()                          // 26
                                            % stats.chrT.m_ind.sn()                          // 27
                                            % (noSP ? "-" : toString(stats.chrT.m_ind.sp())) // 28
                                            % stats.chrT.m_ind.pc()                          // 29
                     ).str());
    o.writer->close();
}

void VDiscover::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);

    o.logInfo("Significance: " + std::to_string(o.sign));
    
    o.logInfo("TP: " + std::to_string(stats.chrT.tps.size()));
    o.logInfo("FP: " + std::to_string(stats.chrT.fps.size()));
    o.logInfo("TN: " + std::to_string(stats.chrT.tns.size()));
    o.logInfo("FN: " + std::to_string(stats.chrT.fns.size()));

    o.info("Generating statistics");

    /*
     * Generating summary statistics
     */
    
    writeSummary("VarDiscover_summary.stats", file, stats, o);

    /*
     * Generating detailed statistics for the sequins
     */
    
    writeQuins("VarDiscover_quins.stats", stats, o);

    /*
     * Generating detailed statistics for the queries
     */
    
    writeQuery("VarDiscover_queries.stats", stats, o);

    switch (o.soft)
    {
        case Software::VarScan:
        {
            /*
             * Generating VarDiscover_ROC.R
             */
            
            o.generate("VarDiscover_ROC.R");
            o.writer->open("VarDiscover_ROC.R");
            o.writer->write(RWriter::createScript("VarDiscover_queries.stats", PlotVROC()));
            o.writer->close();
            
            /*
             * Generating VarDiscover_prob.R
             */
            
            o.generate("VarDiscover_prob.R");
            o.writer->open("VarDiscover_prob.R");
            o.writer->write(RWriter::createScript("VarDiscover_queries.stats", PlotVProb()));
            o.writer->close();

            /*
             * Generating VarDiscover_report.pdf
             */
            
            o.generate("VarDiscover_report.pdf");
            o.report->open("VarDiscover_report.pdf");
            o.report->addTitle("VarDiscover");
            o.report->addFile("VarDiscover_summary.stats");
            o.report->addFile("VarDiscover_quins.stats");
            o.report->addFile("VarDiscover_queries.stats");
            o.report->addFile("VarDiscover_ROC.R");
            o.report->addFile("VarDiscover_Prob.R");

            break;
        }

        default:
        {
            /*
             * Generating a PDF report
             */
            
            o.report->open("VarDiscover_report.pdf");
            o.report->addTitle("VarDiscover");
            o.report->addFile("VarDiscover_summary.stats");
            o.report->addFile("VarDiscover_quins.stats");
            o.report->addFile("VarDiscover_queries.stats");

            break;
        }
    }
}