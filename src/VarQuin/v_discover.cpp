#include "VarQuin/v_discover.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVLOD();

// Defined in resources.cpp
extern Scripts PlotVROC1(); // Somatic (eg: VarScan)

// Defined in resources.cpp
extern Scripts PlotVROC2(); // Germline (eg: GATK)

VDiscover::Stats VDiscover::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VDiscover::Stats stats;
    stats.hist = r.vHist();

    parseVariants(file, o.input, [&](const VariantMatch &m)
    {
        const auto &cID = m.query.cID;

        if (Standard::isSynthetic(cID))
        {
            stats.n_syn++;
            
            /*
             * If no p-value is given (eg: GATK), we'd set it to zero so that the algorithm itself remains unchanged.
             */
            
            const auto p = isnan(m.query.p) ? 0.0 : m.query.p;
            
            // Only matching if the position and alleles agree
            const auto matched = m.match && m.ref && m.alt;
            
            /*
             * Matched by position? reference allele? alternative allele?
             */
            
            if (matched)
            {
                const auto key = var2hash(m.match->id, m.match->type(), m.match->l);
                stats.hist.at(cID).at(key)++;

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
                /*
                 * Variant not found in the reference. This is a FP unless it can be filtered by
                 * the p-value.
                 */

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
        else
        {
            stats.n_gen++;
            stats.geno[cID].vars.push_back(m.query);

            if (m.query.type() == Mutation::SNP)
            {
                stats.geno[cID].snp++;
            }
            else
            {
                stats.geno[cID].ind++;
            }
        }
    });

    /*
     * Sorting out true positives
     */
    
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

    /*
     * Sorting out false positives
     */
    
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

    /*
     * Sorting out true negatives
     */
    
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

    /*
     * Sorting out false negatives
     */
    
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
    stats.chrT.m_snp.nr() = r.countSNPSyn();
    stats.chrT.m_ind.nq() = stats.chrT.dInd();
    stats.chrT.m_ind.nr() = r.countIndSyn();

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
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%";

    o.generate(file);
    o.writer->open("VarDiscover_sequins.csv");
    o.writer->write((boost::format(format) % "ID"
                                           % "Pos"
                                           % "Label"
                                           % "Pval"
                                           % "ReadsR"
                                           % "ReadsV"
                                           % "EFold"
                                           % "EFreq"
                                           % "MFreq"
                                           % "Type").str());

    const auto &r = Standard::instance().r_var;
    
    for (const auto &i : stats.hist)
    {
        const auto cID = i.first;
        
        if (Standard::isSynthetic(cID))
        {
            // For all the variants...
            for (const auto &j : i.second)
            {
                auto key = j.first;
                
                // Detected the variant?
                if (j.second)
                {
                    const auto m = r.findVar(cID, key);
                    assert(m);
                    
                    /*
                     * Now we need to know the label for this reference variant
                     */
                    
                    auto f = [&](const std::vector<VariantMatch> &x, const std::string &label)
                    {
                        for (const auto &k : x)
                        {
                            if (k.match && key == k.match->key())
                            {
                                o.writer->write((boost::format(format) % m->id
                                                                       % m->l.start
                                                                       % label
                                                                       % (isnan(k.query.p) ? "-" : p2str(k.query.p))
                                                                       % k.query.readR
                                                                       % k.query.readV
                                                                       % k.eFold
                                                                       % k.eAllFreq
                                                                       % k.query.alleleFreq()
                                                                       % type2str(m->type())).str());
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
                
                // Failed to detect the variant
                else
                {
                    const auto m = r.findVar(i.first, j.first);
                    assert(m);
                    
                    o.writer->write((boost::format(format) % m->id
                                                           % m->l.start
                                                           % "FN"
                                                           % "NA"
                                                           % "NA"
                                                           % "NA"
                                                           % r.findAFold(m->id)
                                                           % r.findAFreq(m->id)
                                                           % "NA"
                                                           % type2str(m->type())).str());
                }
            }
        }
    }
    
    o.writer->close();
}

static void writeQueries(const FileName &file, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Pos"
                                           % "Label"
                                           % "Ref"
                                           % "Var"
                                           % "EFold"
                                           % "EAllele"
                                           % "Pval"
                                           % "Type").str());

    auto f = [&](const std::vector<VariantMatch> &x, const std::string &label)
    {
        for (const auto &i : x)
        {
            o.writer->write((boost::format(format) % (i.match ? i.match->id : "-")
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

    extern FileName VCFRef();
    extern FileName BedRef();
    
    const auto summary = "-------VarDiscover Output Results\n\n"
                         "-------VarDiscover Output\n\n"
                         "       Reference sequin regions:      %1%\n"
                         "       Reference variant annotations: %2%\n"
                         "       User identified variants:      %3%\n\n"
                         "-------Reference variant annotations\n\n"
                         "       Synthetic: %4% SNPs\n"
                         "       Synthetic: %5% indels\n"
                         "       Synthetic: %6% variants\n\n"
                         "       Genome: %7% SNPs\n"
                         "       Genome: %8% indels\n"
                         "       Genome: %9% variants\n\n"
                         "-------User variant annotations\n\n"
                         "       Synthetic: %10% SNPs\n"
                         "       Synthetic: %11% indels\n"
                         "       Synthetic: %12% variants\n\n"
                         "       Genome: %13% SNPs\n"
                         "       Genome: %14% indels\n"
                         "       Genome: %15% variants\n\n"
                         "-------Identification of synthetic variants\n\n"
                         "       Detected: %16% SNPs\n"
                         "       Detected: %17% indels\n"
                         "       Detected: %18% variants\n\n"
                         "       True Positive:  %19% SNPS\n"
                         "       True Positive:  %20% indels\n"
                         "       True Positive:  %21% variants\n\n"
                         "       False Positive: %22% SNPs\n"
                         "       False Positive: %23% indels\n"
                         "       False Positive: %24% variants\n\n"
                         "       False Negative: %25% SNPs\n"
                         "       False Negative: %26% indels\n"
                         "       False Negative: %27% variants\n\n"
                         "-------Diagnostic Performance\n\n"
                         "       *Variants\n"
                         "       Sensitivity: %28%\n"
                         "       Specificity: %29%\n"
                         "       Precision:   %30%\n"
                         "       FDR Rate:    %31%\n\n"
                         "       *SNVs\n"
                         "       Sensitivity: %32%\n"
                         "       Specificity: %33%\n"
                         "       Precision:   %34%\n"
                         "       FDR Rate:    %35%\n\n"
                         "       *Indels\n"
                         "       Sensitivity: %36%\n"
                         "       Specificity: %37%\n"
                         "       Precision:   %38%\n"
                         "       FDR Rate:    %39%\n";
    o.generate(file);
    o.writer->open("VarDiscover_summary.stats");
    o.writer->write((boost::format(summary) % BedRef()
                                            % VCFRef()
                                            % src
                                            % r.countSNPSyn()
                                            % r.countIndSyn()
                                            % (r.countSNPSyn() + r.countIndSyn())
                                            % r.countSNPGen()
                                            % r.countIndGen()
                                            % (r.countSNPGen() + r.countIndGen())
                                            % stats.chrT.dSNP()
                                            % stats.chrT.dInd()
                                            % stats.chrT.dTot()
                                            % stats.countSNPGeno()
                                            % stats.countIndGeno()
                                            % stats.countVarGeno()
                                            % stats.chrT.dSNP()
                                            % stats.chrT.dInd()
                                            % stats.chrT.dTot()
                                            % stats.chrT.tpSNP()
                                            % stats.chrT.tpInd()
                                            % stats.chrT.tpTot()
                                            % stats.chrT.fpSNP()
                                            % stats.chrT.fpInd()
                                            % stats.chrT.fpTot()
                                            % stats.chrT.fnSNP()
                                            % stats.chrT.fnInd()     // 26
                                            % stats.chrT.fnTot()     // 27
                                            % stats.chrT.m.sn()      // 28
                                            % stats.chrT.m.sp()      // 29
                                            % stats.chrT.m.pc()      // 30
                                            % stats.chrT.m.fdr()     // 31
                                            % stats.chrT.m_snp.sn()  // 32
                                            % stats.chrT.m_snp.sp()  // 33
                                            % stats.chrT.m_snp.pc()  // 34
                                            % stats.chrT.m_snp.fdr() // 35
                                            % stats.chrT.m_ind.sn()  // 36
                                            % stats.chrT.m_ind.sp()  // 37
                                            % stats.chrT.m_ind.pc()  // 38
                                            % stats.chrT.m_snp.fdr() // 39
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
     * Generating VarDiscover_quins.stats
     */
    
    writeQuins("VarDiscover_sequins.csv", stats, o);

    /*
     * Generating VarDiscover_summary.stats
     */
    
    writeSummary("VarDiscover_summary.stats", file, stats, o);
    
    /*
     * Generating VarDiscover_queries.stats
     */
    
    writeQueries("VarDiscover_queries.stats", stats, o);
    
    /*
     * Generating VarDiscover_ROC.R
     */
    
    o.generate("VarDiscover_ROC.R");
    o.writer->open("VarDiscover_ROC.R");
    o.writer->write(RWriter::createScript("VarDiscover_queries.stats", PlotVROC1()));
    o.writer->close();
    
    /*
     * Generating VarDiscover_LOD.R
     */
    
    o.generate("VarDiscover_LOD.R");
    o.writer->open("VarDiscover_LOD.R");
    o.writer->write(RWriter::createScript("VarDiscover_queries.stats", PlotVLOD()));
    o.writer->close();
    
    /*
     * Generating VarDiscover_report.pdf
     */
    
    o.report->open("VarDiscover_report.pdf");
    o.report->addTitle("VarDiscover");
    o.report->addFile("VarDiscover_summary.stats");
    o.report->addFile("VarDiscover_sequins.csv");
    o.report->addFile("VarDiscover_queries.stats");
    o.report->addFile("VarDiscover_ROC.R");
    o.report->addFile("VarDiscover_LOD.R");
}