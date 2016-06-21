#include "VarQuin/v_discover.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVLOD();

// Defined in resources.cpp
extern Scripts PlotSomaticROC();

// Defined in resources.cpp
extern Scripts PlotGermlineROC();

// Defined in standard.cpp
extern bool IsFlatMix();

VDiscover::Stats VDiscover::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VDiscover::Stats stats;
    stats.hist = r.vHist();

    for (const auto &i : stats.hist)
    {
        stats.data[i.first];
    }

    o.info("Reading VCF inputs");
    stats.vData = vcfData(file, o.input);

    o.info("Parsing VCF inputs");
    parseVariants(file, o.input, [&](const VariantMatch &m)
    {
        const auto &cID = m.query.cID;

        auto f = [&]()
        {
            /*
             * If no p-value is given (eg: GATK), we'd set it to zero so that the algorithm itself
             * remains unchanged.
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
                    stats.data[cID].tps.push_back(m);
                    stats.data[cID].tps_[key] = &stats.data[cID].tps.back();
                }
                else
                {
                    stats.data[cID].fns.push_back(m);
                    stats.data[cID].fns_[key] = &stats.data[cID].fns.back();
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
                    //stats.data[cID].fps_.insert(key);
                    stats.data[cID].fps.push_back(m);
                }
                else
                {
                    //stats.data[cID].tns_.insert(key);
                    stats.data[cID].tns.push_back(m);
                }
            }
        };

        if (Standard::isSynthetic(cID))
        {
            stats.n_syn++;
            f();
        }
        else
        {
            stats.n_gen++;
            
            if (Standard::isGenomic(cID))
            {
                f();
            }
        }
    });

    for (const auto &i : stats.data)
    {
        auto &x = stats.data[i.first];
        
        for (const auto &j : i.second.tps)
        {
            i.second.m.tp()++;
            
            switch (j.query.type())
            {
                case Mutation::SNP:       { x.m_snp.tp()++; break; }
                case Mutation::Deletion:
                case Mutation::Insertion: { x.m_ind.tp()++; break; }
            }
        }
        
        for (const auto &j : i.second.fps)
        {
            i.second.m.fp()++;
            
            switch (j.query.type())
            {
                case Mutation::SNP:       { x.m_snp.fp()++; break; }
                case Mutation::Deletion:
                case Mutation::Insertion: { x.m_ind.fp()++; break; }
            }
        }
        
        for (const auto &j : i.second.tns)
        {
            i.second.m.tn()++;
            
            switch (j.query.type())
            {
                case Mutation::SNP:       { x.m_snp.tn()++; break; }
                case Mutation::Deletion:
                case Mutation::Insertion: { x.m_ind.tn()++; break; }
            }
        }
        
        for (const auto &j : i.second.fns)
        {
            i.second.m.fn()++;
            
            switch (j.query.type())
            {
                case Mutation::SNP:       { x.m_snp.fn()++; break; }
                case Mutation::Deletion:
                case Mutation::Insertion: { x.m_ind.fn()++; break; }
            }
        }

        x.m_snp.nq() = x.dSNP();
        x.m_snp.nr() = r.countSNP(i.first);
        x.m_ind.nq() = x.dInd();
        x.m_ind.nr() = r.countInd(i.first);

        x.m.nq() = x.m_snp.nq() + x.m_ind.nq();
        x.m.nr() = x.m_snp.nr() + x.m_ind.nr();
    }
    
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
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";

    o.generate(file);
    o.writer->open("VarDiscover_sequins.csv");
    o.writer->write((boost::format(format) % "ID"
                                           % "Pos"
                                           % "Label"
                                           % "Pval"
                                           % "ReadR"
                                           % "ReadV"
                                           % "Depth"
                                           % "EFold"
                                           % "EFreq"
                                           % "MFreq"
                                           % "Type").str());

    const auto &r = Standard::instance().r_var;
    
    for (const auto &i : stats.hist)
    {
        if (!Standard::isSynthetic(i.first))
        {
            continue;
        }
        
        std::cout << i.second.size() << std::endl;
        
//        for (const auto &j : i.second)
        {
            const auto cID = i.first;
            
//            if (Standard::isSynthetic(cID))
            {
                int p = 0;
                
                // Search all query variants...
                for (const auto &j : i.second)
                {
                    std::cout << p++ << std::endl;
                    
                    auto key = j.first;
                    
                    // Detected the sequin?
                    if (j.second)
                    {
                        const auto m = r.findVar(cID, key);
                        assert(m);
                        
                        /*
                         * Now we need to know the label for this reference variant
                         */
                        
                        auto f = [&](const std::map<long, VariantMatch *> &x, const std::string &label)
                        {
//                            for (const auto &k : x)
                            {
                                if (x.count(key)/* k.match && key == k.match->key() */)
                                {
                                    const auto t = x.at(key);
                                    
                                    o.writer->write((boost::format(format)
                                                     % m->id
                                                     % m->l.start
                                                     % label
                                                     % (isnan(t->query.p) ? "-" : p2str(t->query.p))
                                                     % t->query.readR
                                                     % t->query.readV
                                                     % t->query.cov
                                                     % t->eFold
                                                     % t->eAllFreq
                                                     % t->query.alleleFreq()
                                                     % type2str(m->type())).str());
                                    return true;
                                }
                            }
                            
                            return false;
                        };
                        
                        if (!f(stats.data.at(ChrT).tps_, "TP") &&
                            //!f(stats.data.at(ChrT).fps, "FP") &&
                            //!f(stats.data.at(ChrT).tns, "TN") &&
                            !f(stats.data.at(ChrT).fns_, "FN")
                            )
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
                                         % "NA"
                                         % r.findAFold(m->id)
                                         % r.findAFreq(m->id)
                                         % "NA"
                                         % type2str(m->type())).str());
                    }
                }
            }
        }
    }
    
    o.writer->close();
}

static void writeQueries(const FileName &file, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Pos"
                                           % "Label"
                                           % "Ref"
                                           % "Var"
                                           % "Depth"
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
                                                   % i.query.cov
                                                   % i.eFold
                                                   % i.eAllFreq
                                                   % (isnan(i.query.p) ? "-" : p2str(i.query.p))
                                                   % type2str(i.query.type())).str());
        }
    };

    f(stats.data.at(ChrT).tps, "TP");
    f(stats.data.at(ChrT).fps, "FP");
    f(stats.data.at(ChrT).tns, "TN");
    f(stats.data.at(ChrT).fns, "FN");

    o.writer->close();
}

static void writeSummary(const FileName &file, const FileName &src, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;

    extern FileName VCFRef();
    extern FileName BedRef();
    
    const auto summary = "-------VarDiscover Output Results\n\n"
                         "-------VarDiscover Output\n\n"
                         "       Reference variant annotations: %1%\n"
                         "       User identified variants:      %2%\n\n"
                         "-------Reference variant annotations\n\n"
                         "       Synthetic: %3% SNPs\n"
                         "       Synthetic: %4% indels\n"
                         "       Synthetic: %5% variants\n\n"
                         "       Genome: %6% SNPs\n"
                         "       Genome: %7% indels\n"
                         "       Genome: %8% variants\n\n"
                         "-------User identified variants\n\n"
                         "       Synthetic: %9% SNPs\n"
                         "       Synthetic: %10% indels\n"
                         "       Synthetic: %11% variants\n\n"
                         "       Genome: %12% SNPs\n"
                         "       Genome: %13% indels\n"
                         "       Genome: %14% variants\n\n"
                         "-------Identification of synthetic variants\n\n"
                         "       True Positive:  %15% SNPS\n"
                         "       True Positive:  %16% indels\n"
                         "       True Positive:  %17% variants\n\n"
                         "       False Positive: %18% SNPs\n"
                         "       False Positive: %19% indels\n"
                         "       False Positive: %20% variants\n\n"
                         "       False Negative: %21% SNPs\n"
                         "       False Negative: %22% indels\n"
                         "       False Negative: %23% variants\n\n"
                         "-------Diagnostic Performance (Synthetc)\n\n"
                         "       *Variants\n"
                         "       Sensitivity: %24%\n"
                         "       Specificity: %25%\n"
                         "       Precision:   %26%\n"
                         "       FDR Rate:    %27%\n\n"
                         "       *SNVs\n"
                         "       Sensitivity: %28%\n"
                         "       Specificity: %29%\n"
                         "       Precision:   %30%\n"
                         "       FDR Rate:    %31%\n\n"
                         "       *Indels\n"
                         "       Sensitivity: %32%\n"
                         "       Specificity: %33%\n"
                         "       Precision:   %34%\n"
                         "       FDR Rate:    %35%\n\n"
                         "-------Identification of genomic variants\n\n"
                         "       True Positive:  %36% SNPS\n"
                         "       True Positive:  %37% indels\n"
                         "       True Positive:  %38% variants\n\n"
                         "       False Positive: %39% SNPs\n"
                         "       False Positive: %40% indels\n"
                         "       False Positive: %41% variants\n\n"
                         "       False Negative: %42% SNPs\n"
                         "       False Negative: %43% indels\n"
                         "       False Negative: %44% variants\n\n"
                         "-------Diagnostic Performance (Genome)\n\n"
                         "       *Variants\n"
                         "       Sensitivity: %45%\n"
                         "       Specificity: %46%\n"
                         "       Precision:   %47%\n"
                         "       FDR Rate:    %48%\n\n"
                         "       *SNVs\n"
                         "       Sensitivity: %49%\n"
                         "       Specificity: %50%\n"
                         "       Precision:   %51%\n"
                         "       FDR Rate:    %52%\n\n"
                         "       *Indels\n"
                         "       Sensitivity: %53%\n"
                         "       Specificity: %54%\n"
                         "       Precision:   %55%\n"
                         "       FDR Rate:    %56%\n";
    o.generate(file);
    o.writer->open("VarDiscover_summary.stats");
    o.writer->write((boost::format(summary) % VCFRef()
                                            % src
                                            % r.countSNPSyn()
                                            % r.countIndSyn()
                                            % (r.countSNPSyn() + r.countIndSyn())
                                            % r.countSNPGen()
                                            % r.countIndGen()
                                            % (r.countSNPGen() + r.countIndGen())

                                            % stats.vData.countSNPSyn()
                                            % stats.vData.countIndSyn() // 10
                                            % stats.vData.countVarSyn() // 11
                                            % stats.vData.countSNPGen() // 12
                                            % stats.vData.countIndGen() // 13
                                            % stats.vData.countVarGen() // 14
                                            % stats.countSNP_TP_Syn()     // 15
                                            % stats.countInd_TP_Syn()     // 16
                                            % stats.countVar_TP_Syn()     // 17
                                            % stats.countSNP_FP_Syn() // 18
                                            % stats.countInd_FP_Syn() // 19
                                            % stats.countVar_FP_Syn() // 20
                                            % stats.countSNP_FN_Syn()     // 21
                                            % stats.countInd_FN_Syn()     // 22
                                            % stats.countVar_FN_Syn()     // 23
                     // 0.893773
                                            % stats.countVarSN_Syn()//stats.data.at(ChrT).m.sn() //      // 24
                                            % stats.countVarSP_Syn()      // 25
                                            % stats.countVarPC_Syn()      // 26
                                            % (1-stats.countVarPC_Syn())     // 27
                                            % stats.countSNPSN_Syn()  // 28
                                            % stats.countSNPSP_Syn()  // 29
                                            % stats.countSNPPC_Syn()  // 30
                                            % (1-stats.countSNPPC_Syn()) // 31
                                            % stats.countIndSN_Syn()  // 32
                                            % stats.countIndSP_Syn()  // 33
                                            % stats.countIndPC_Syn()  // 34
                                            % (1-stats.countIndPC_Syn()) // 35

                     % stats.countSNP_TP_Gen()     // 36
                     % stats.countInd_TP_Gen()     // 37
                     % stats.countVar_TP_Gen()     // 38
                     % stats.countSNP_FP_Gen() // 39
                     % stats.countInd_FP_Gen() // 40
                     % stats.countVar_FP_Gen() // 41
                     % stats.countSNP_FN_Gen()     // 42
                     % stats.countInd_FN_Gen()     // 43
                     % stats.countVar_FN_Gen()     // 44
                    
                     % stats.countVarSN_Gen() // 45
                     % stats.countVarSP_Gen() // 46
                     % stats.countVarPC_Gen() // 47
                     % (1-stats.countVarPC_Gen()) // 48
                     % stats.countSNPSN_Gen() // 49
                     % stats.countSNPSP_Gen() // 50
                     % stats.countSNPPC_Gen() // 51
                     % (1-stats.countSNPPC_Gen()) // 52
                     % stats.countIndSN_Gen() // 53
                     % stats.countIndSP_Gen() // 54
                     % stats.countIndPC_Gen() // 55
                     % (1-stats.countIndPC_Gen()) // 56
                     
                     
                     ).str());
    o.writer->close();
}

void VDiscover::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);

    o.logInfo("Significance: " + std::to_string(o.sign));
    
    o.logInfo("TP: " + std::to_string(stats.data.at(ChrT).tps.size()));
    o.logInfo("FP: " + std::to_string(stats.data.at(ChrT).fps.size()));
    o.logInfo("TN: " + std::to_string(stats.data.at(ChrT).tns.size()));
    o.logInfo("FN: " + std::to_string(stats.data.at(ChrT).fns.size()));

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
    
    if (IsFlatMix())
    {
        o.writer->write(RWriter::createScript("VarDiscover_queries.stats", PlotGermlineROC()));
    }
    else
    {
        o.writer->write(RWriter::createScript("VarDiscover_queries.stats", PlotSomaticROC()));
    }
    
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