#include "VarQuin/v_freq.hpp"
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
            
            //const auto p = 0; //isnan(m.query.p) ? 0.0 : m.query.p;
            
            // Only matching if the position and alleles agree
            const auto matched = m.match && m.ref && m.alt;
            
            /*
             * Matched by position? reference allele? alternative allele?
             */
            
            if (matched)
            {
                const auto key = var2hash(m.match->id, m.match->type(), m.match->l);
                stats.hist.at(cID).at(key)++;
                
               // if (p <= o.sign)
                {
                    stats.data[cID].tps.push_back(m);
                    stats.data[cID].tps_[key] = stats.data[cID].tps.back();
                }
//                else
//                {
//                    throw "Not Implemented";
//                    //stats.data[cID].fns.push_back(m);
//                    //stats.data[cID].fns_[key] = &stats.data[cID].fns.back();
//                }
            }
            else
            {
                /*
                 * Variant not found in the reference. This is a FP unless it can be filtered by
                 * the p-value.
                 */
                
                //if (p <= o.sign)
                {
                    //stats.data[cID].fps_.insert(key);
                    stats.data[cID].fps.push_back(m);
                }
//                else
//                {
//                    throw "?????";
//                    //stats.data[cID].tns_.insert(key);
//                    //stats.data[cID].tns.push_back(m);
//                }
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

    o.info("Collecting statistics");
    
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
        
        x.m_snp.nq() = x.dSNP();
        x.m_snp.nr() = r.countSNP(i.first);
        x.m_ind.nq() = x.dInd();
        x.m_ind.nr() = r.countInd(i.first);

        x.m.nq() = x.m_snp.nq() + x.m_ind.nq();
        x.m.nr() = x.m_snp.nr() + x.m_ind.nr();

        x.m.fn()     = x.m.nr()     - x.m.tp();
        x.m_snp.fn() = x.m_snp.nr() - x.m_snp.tp();
        x.m_ind.fn() = x.m_ind.nr() - x.m_ind.tp();

        assert(x.m.nr() >= x.m.fn());
        assert(x.m_snp.nr() >= x.m_snp.fn());
        assert(x.m_ind.nr() >= x.m_ind.fn());
    }
    
    return stats;
}

static void writeQuins(const FileName &file,
                       const VDiscover::Stats &stats,
                       const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;

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
    for (const auto &i : stats.hist)
    {
        if (!Standard::isSynthetic(i.first))
        {
            continue;
        }
        
        const auto cID = i.first;
        
        // Search all query variants...
        for (const auto &j : i.second)
        {
            auto key = j.first;
            
            // Detected the sequin?
            if (j.second)
            {
                const auto m = r.findVar(cID, key);
                assert(m);
                
                /*
                 * Now we need to know the label for this reference variant
                 */
                
                auto f = [&](const std::map<long, VariantMatch> &x, const std::string &label)
                {
                    if (x.count(key))
                    {
                        const auto &t = x.at(key);
                        
                        o.writer->write((boost::format(format) % m->id
                                                               % m->l.start
                                                               % label
                                                               % (isnan(t.query.p) ? "-" : p2str(t.query.p))
                                                               % t.query.readR
                                                               % t.query.readV
                                                               % t.query.depth
                                                               % t.eFold
                                                               % t.eAllFreq
                                                               % t.query.alleleFreq()
                                                               % type2str(m->type())).str());
                        return true;
                    }
                    
                    return false;
                };

                if (!f(stats.data.at(i.first).tps_, "TP") && !f(stats.data.at(i.first).fns_, "FN"))
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
            const auto eFold = (label == "FP" ? NAN : i.eFold);
            const auto eAlFq = (label == "FP" ? NAN : i.eAllFreq);
            
            o.writer->write((boost::format(format) % (i.match ? i.match->id : "-")
                                                   % i.query.l.start
                                                   % label
                                                   % i.query.readR
                                                   % i.query.readV
                                                   % i.query.depth
                                                   % eFold
                                                   % eAlFq
                                                   % (isnan(i.query.p) ? "-" : p2str(i.query.p))
                                                   % type2str(i.query.type())).str());
        }
    };

    for (const auto &i : stats.data)
    {
        if (Standard::isSynthetic(i.first))
        {
            f(i.second.tps, "TP");
            f(i.second.fps, "FP");
        }
    }

    o.writer->close();
}

static void writeSummary(const FileName &file, const FileName &src, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;

    extern FileName VCFRef();
    extern FileName BedRef();
    
    const auto hasGen = stats.data.size() > 1;
    
    #define G(x) (hasGen ? "-" : toString(x))
    
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
                         "-------Diagnostic Performance (Synthetic)\n\n"
                         "       *Variants\n"
                         "       Sensitivity: %24$.4f\n"
                         "       Precision:   %25$.4f\n"
                         "       FDR Rate:    %26$.4f\n\n"
                         "       *SNVs\n"
                         "       Sensitivity: %27$.4f\n"
                         "       Precision:   %28$.4f\n"
                         "       FDR Rate:    %29$.4f\n\n"
                         "       *Indels\n"
                         "       Sensitivity: %30$.4f\n"
                         "       Precision:   %31$.4f\n"
                         "       FDR Rate:    %32$.4f\n\n"
                         "-------Identification of genomic variants\n\n"
                         "       True Positive:  %33% SNPS\n"
                         "       True Positive:  %34% indels\n"
                         "       True Positive:  %35% variants\n";
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
                                            % stats.vData.countSNPSyn()  // 9
                                            % stats.vData.countIndSyn()  // 10
                                            % stats.vData.countVarSyn()  // 11
                                            % stats.vData.countSNPGen()  // 12
                                            % stats.vData.countIndGen()  // 13
                                            % stats.vData.countVarGen()  // 14
                                            % stats.countSNP_TP_Syn()    // 15
                                            % stats.countInd_TP_Syn()    // 16
                                            % stats.countVar_TP_Syn()    // 17
                                            % stats.countSNP_FP_Syn()    // 18
                                            % stats.countInd_FP_Syn()    // 19
                                            % stats.countVar_FP_Syn()    // 20
                                            % stats.countSNP_FN_Syn()    // 21
                                            % stats.countInd_FN_Syn()    // 22
                                            % stats.countVar_FN_Syn()    // 23
                                            % stats.countVarSN_Syn()     // 24
                                            % stats.countVarPC_Syn()     // 25
                                            % (1-stats.countVarPC_Syn()) // 26
                                            % stats.countSNPSN_Syn()     // 27
                                            % stats.countSNPPC_Syn()     // 28
                                            % (1-stats.countSNPPC_Syn()) // 29
                                            % stats.countIndSN_Syn()     // 30
                                            % stats.countIndPC_Syn()     // 31
                                            % (1-stats.countIndPC_Syn()) // 32
                                            % G(stats.countSNP_TP_Gen()) // 33
                                            % G(stats.countInd_TP_Gen()) // 34
                                            % G(stats.countVar_TP_Gen()) // 35
                     ).str());
    o.writer->close();
}

void VDiscover::report(const FileName &file, const Options &o)
{
    // Statistics for the variants
    const auto stats = analyze(file, o);
    
    o.info("TP: " + std::to_string(stats.countVar_TP_Syn()));
    o.info("FP: " + std::to_string(stats.countVar_FP_Syn()));
    o.info("FN: " + std::to_string(stats.countVar_FN_Syn()));

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