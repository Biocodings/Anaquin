#include "VarQuin/v_freq.hpp"
#include "VarQuin/v_discover.hpp"
#include <fstream>
using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVLOD();

// Defined in resources.cpp
extern Scripts PlotVROC();

// Defined in resources.cpp
extern Scripts PlotGermlineROC();

static Counts __countD__ = 0;
static Counts __countP__ = 0;

std::ofstream out("fuck.txt");

struct VDiscoverImpl : public VCFDataUser
{
    VDiscover::Stats *stats;
    
    void variantProcessed(const ParserVCF::Data &x, const ParserProgress &p) override
    {
        const auto &r = Standard::instance().r_var;
        
        VariantMatch m;

        auto match = [&](const CalledVariant &query)
        {
            m.query = query;
            m.match = nullptr;
            
            const auto isSyn = Standard::isSynthetic(query.cID);
            
            if (isSyn || Standard::isGenomic(query.cID))
            {
                // Can we match by position?
                m.match = r.findVar(query.cID, query.l);
                
                if (m.match)
                {
                    m.ref = m.match->ref == query.ref;
                    m.alt = m.match->alt == query.alt;
                }
                
                const auto perfect = m.match && m.ref && m.alt;
                
                if (isSyn && perfect)
                {
                    m.eFold    = r.findAFold(baseID(m.match->id));
                    m.eAllFreq = r.findAFreq(baseID(m.match->id));
                }
                else
                {
                    m.eFold    = NAN;
                    m.eAllFreq = NAN;
                }
            }
            
            return m.match;
        };
     
        match(x);
        
        const auto &cID = m.query.cID;
        
        auto f = [&]()
        {
            if (!isnan(m.query.p))     { __countP__++; }
            if (!isnan(m.query.depth)) { __countD__++; }

            /*
             * If no p-value is given (eg: GATK), we'd set it to zero so that the algorithm itself
             * remains unchanged.
             */
            
            // Only matching if the position and alleles agree
            const auto matched = m.match && m.ref && m.alt;
            
            /*
             * Matched by position? reference allele? alternative allele?
             */
            
            if (matched)
            {
                const auto key = m.match->key();
                stats->hist.at(cID).at(key)++;
                
                stats->data[cID].tps.push_back(m);
                stats->data[cID].tps_[key] = stats->data[cID].tps.back();
                
                //stats.hist.at(cID).at(key)++;
                //stats.hist.at(m.match->id)++;
                
                stats->data.at(cID).af = m.query.alleleFreq();
                
                if (Standard::isSynthetic(cID))
                {
                    const auto exp = r.findAFreq(baseID(m.match->id));
                    const auto obs = m.query.alleleFreq();

                    out << m.match->l.start << "\t" << exp << "\t" << obs << std::endl;
                    
                    
                    // Eg: 2821292107
                    const auto id = toString(key);
                    
                    // Add for all variants
                    stats->vars.add(id, exp, obs);
                    
                    switch (m.query.type())
                    {
                        case Mutation::SNP:       { stats->snp.add(id, exp, obs); break; }
                        case Mutation::Deletion:
                        case Mutation::Insertion: { stats->ind.add(id, exp, obs); break; }
                    }
                    
                    stats->readR[key] = m.query.readR;
                    stats->readV[key] = m.query.readV;
                    stats->depth[key] = m.query.depth;
                    
                    if (isnan(stats->vars.limit.abund) || exp < stats->vars.limit.abund)
                    {
                        stats->vars.limit.id = m.match->id;
                        stats->vars.limit.abund = exp;
                    }
                }
            }
            else
            {
                /*
                 * Variant not found in the reference. This is a FP unless it can be filtered by
                 * the p-value.
                 */
                
                //stats.data[cID].fps_.insert(key);
                stats->data[cID].fps.push_back(m);
            }
        };
        
        if (Standard::isSynthetic(cID))
        {
            stats->n_syn++;
            f();
        }
        else
        {
            stats->n_gen++;
            
            if (Standard::isGenomic(cID))
            {
                f();
            }
        }
    }
};

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

    VDiscoverImpl impl;
    impl.stats = &stats;

    // Read the input variant file and process them
    stats.vData = vcfData(file, o.input, &impl);

    o.info("Aggregating statistics");
    out.close();
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
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%";
    
    o.generate(file);
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

    for (const auto &i : stats.data)
    {
        if (Standard::isSynthetic(i.first))
        {
            const auto &cID = i.first;
            
            // We'll need it to search for the sequin where the FPs are
            const auto inters = r.mInters(cID);
            
            assert(inters.size());

            auto f = [&](const std::vector<VariantMatch> &x, const std::string &label)
            {
                for (const auto &i : x)
                {
                    auto sID = (i.match ? i.match->id : "-");
                    
                    auto eFold = i.eFold;
                    auto eAllF = i.eAllFreq;
                    
                    if (label == "FP")
                    {
                        assert(isnan(eFold));
                        assert(isnan(eAllF));
                        
                        const auto m = inters.contains(i.query.l);
                        
                        // Can we find the corresponding region for the FP?
                        if (m)
                        {
                            sID = m->id();

                            // It has to be sequin name (eg: D_3_12)
                            assert(!sID.empty());
                            
                           eFold = r.findAFold(sID);
                           eAllF = r.findAFreq(sID);
                        }
                    }
                    
                    const auto pval  = (isnan(i.query.p) ? "-" : p2str(i.query.p));
                    
                    o.writer->write((boost::format(format) % sID
                                                           % i.query.l.start
                                                           % label
                                                           % i.query.readR
                                                           % i.query.readV
                                                           % i.query.depth
                                                           % eFold
                                                           % eAllF
                                                           % pval
                                                           % type2str(i.query.type())).str());
                }
            };
            
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
    extern FileName MixRef();

    std::cout << stats.vars.size() << std::endl;
    
    const auto lm = stats.vars.linear(true);
    
    // Calcluate the quantification point with logarithm
    auto ms = stats.vars.limitQuant(true);
    
    // Remember the break-point is on the log2-scale, we'll need to convert it back
    ms.b = pow(2, ms.b);
    
    Counts n_below = 0;
    Counts n_above = 0;
    
    for (const auto &i : stats.data)
    {
        if (!Standard::isSynthetic(i.first))
        {
            if (i.second.af >= ms.b)
            {
                n_above++;
            }
            else
            {
                n_below++;
            }
        }
    }
    
    const auto mm = r.findVar(ChrT, stol(ms.id));
    assert(mm);
    
    const auto summary = "-------VarDiscover Output Results\n\n"
                         "-------VarDiscover Output\n\n"
                         "       Reference variant annotations:   %1%\n"
                         "       Referene coordinate annotations: %2%\n"
                         "       Sequin mixture file:             %3%\n"
                         "       User identified variants:        %4%\n\n"
                         "-------Reference variant annotations\n\n"
                         "       Synthetic: %5% SNPs\n"
                         "       Synthetic: %6% indels\n"
                         "       Synthetic: %7% variants\n\n"
                         "       Genome:    %8% SNPs\n"
                         "       Genome:    %9% indels\n"
                         "       Genome:    %10% variants\n\n"
                         "-------User identified variants\n\n"
                         "       Synthetic: %11% SNPs\n"
                         "       Synthetic: %12% indels\n"
                         "       Synthetic: %13% variants\n\n"
                         "       Detection Sensitivity: %14% (attomol/ul) (%15%)\n\n"
                         "       Genome: %16% SNPs\n"
                         "       Genome: %17% indels\n"
                         "       Genome: %18% variants\n\n"
                         "-------Identification of synthetic variants\n\n"
                         "       True Positive:  %19% SNPS\n"
                         "       True Positive:  %20% indels\n"
                         "       True Positive:  %21% variants\n\n"
                         "       False Positive: %22% SNPs\n"
                         "       False Positive: %23% indels\n"
                         "       False Positive: %24% variants\n\n"
                         "       False Negative: %25% SNPs\n"
                         "       False Negative: %26% indels\n"
                         "       False Negative: %27% variants\n\n"
                         "-------Diagnostic Performance (Synthetic)\n\n"
                         "       *Variants\n"
                         "       Sensitivity: %28$.4f\n"
                         "       Precision:   %29$.4f\n"
                         "       FDR Rate:    %30$.4f\n\n"
                         "       *SNVs\n"
                         "       Sensitivity: %31$.4f\n"
                         "       Precision:   %32$.4f\n"
                         "       FDR Rate:    %33$.4f\n\n"
                         "       *Indels\n"
                         "       Sensitivity: %34$.4f\n"
                         "       Precision:   %35$.4f\n"
                         "       FDR Rate:    %36$.4f\n\n"
                         "-------Limit of Quantification (LOQ)\n"
                         "      *Estimated by piecewise segmented regression\n\n"
                         "       Break: %37% attomol/ul (%38%)\n\n"
                         "      *Below LOQ\n"
                         "       Intercept:   %39%\n"
                         "       Slope:       %40%\n"
                         "       Correlation: %41%\n"
                         "       R2:          %42%\n"
                         "       Genome:      %43%\n\n"
                         "      *Above LOQ\n"
                         "       Intercept:   %44%\n"
                         "       Slope:       %45%\n"
                         "       Correlation: %46%\n"
                         "       R2:          %47%\n"
                         "       Genome:      %48%\n\n"
                         "-------Overall linear regression (log2 scale)\n\n"
                         "      Correlation: %49%\n"
                         "      Slope:       %50%\n"
                         "      R2:          %51%\n"
                         "      F-statistic: %52%\n"
                         "      P-value:     %53%\n"
                         "      SSM:         %54%, DF: %55%\n"
                         "      SSE:         %56%, DF: %57%\n"
                         "      SST:         %58%, DF: %59%\n";
    o.generate(file);
    o.writer->open("VarDiscover_summary.stats");
    o.writer->write((boost::format(summary) % VCFRef()                   // 1
                                            % BedRef()                   // 2
                                            % MixRef()                   // 3
                                            % src                        // 4
                                            % r.countSNPSyn()            // 5
                                            % r.countIndSyn()            // 6
                                            % (r.countSNPSyn() + r.countIndSyn())
                                            % r.countSNPGen()            // 8
                                            % r.countIndGen()            // 9
                                            % (r.countSNPGen() + r.countIndGen())
                                            % stats.vData.countSNPSyn()  // 11
                                            % stats.vData.countIndSyn()  // 12
                                            % stats.vData.countVarSyn()  // 13
                                            % stats.vars.limit.abund     // 14
                                            % stats.vars.limit.id        // 15
                                            % stats.vData.countSNPGen()  // 16
                                            % stats.vData.countIndGen()  // 17
                                            % stats.vData.countVarGen()  // 18
                                            % stats.countSNP_TP_Syn()    // 19
                                            % stats.countInd_TP_Syn()    // 20
                                            % stats.countVar_TP_Syn()    // 21
                                            % stats.countSNP_FP_Syn()    // 22
                                            % stats.countInd_FP_Syn()    // 23
                                            % stats.countVar_FP_Syn()    // 24
                                            % stats.countSNP_FN_Syn()    // 25
                                            % stats.countInd_FN_Syn()    // 26
                                            % stats.countVar_FN_Syn()    // 27
                                            % stats.countVarSN_Syn()     // 28
                                            % stats.countVarPC_Syn()     // 29
                                            % (1-stats.countVarPC_Syn()) // 30
                                            % stats.countSNPSN_Syn()     // 31
                                            % stats.countSNPPC_Syn()     // 32
                                            % (1-stats.countSNPPC_Syn()) // 33
                                            % stats.countIndSN_Syn()     // 34
                                            % stats.countIndPC_Syn()     // 35
                                            % (1-stats.countIndPC_Syn()) // 36
                                            % ms.b                       // 37
                                            % mm->id                     // 38
                                            % ms.lInt                    // 39
                                            % ms.lSl                     // 40
                                            % ms.lr                      // 41
                                            % ms.lR2                     // 42
                                            % n_above                    // 43
                                            % ms.rInt                    // 44
                                            % ms.rSl                     // 45
                                            % ms.rr                      // 46
                                            % ms.rR2                     // 47
                                            % n_below                    // 48
                                            % lm.r                       // 49
                                            % lm.m                       // 50
                                            % lm.R2                      // 51
                                            % lm.F                       // 52
                                            % lm.p                       // 53
                                            % lm.SSM                     // 54
                                            % lm.SSM_D                   // 55
                                            % lm.SSE                     // 56
                                            % lm.SSE_D                   // 57
                                            % lm.SST                     // 58
                                            % lm.SST_D                   // 59
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
     * Generating VarDiscover_queries.csv
     */
    
    writeQueries("VarDiscover_queries.csv", stats, o);
    
    /*
     * Generating VarDiscover_ROC.R
     */
    
    o.generate("VarDiscover_ROC.R");
    o.writer->open("VarDiscover_ROC.R");
    
    if (__countP__ >= __countD__)
    {
        o.info("P-value for scoring");
        o.writer->write(RWriter::createVROC("VarDiscover_queries.csv", "1-data$Pval"));
    }
    else
    {
        o.info("Depth for scoring");
        o.writer->write(RWriter::createVROC("VarDiscover_queries.csv", "data$Depth"));
    }
    
    o.writer->close();
    
    /*
     * Generating VarDiscover_LOD.R
     */
    
    o.generate("VarDiscover_LOD.R");
    o.writer->open("VarDiscover_LOD.R");
    o.writer->write(RWriter::createScript("VarDiscover_queries.csv", PlotVLOD()));
    o.writer->close();
    
    /*
     * Generating VarDiscover_report.pdf
     */
    
    o.report->open("VarDiscover_report.pdf");
    o.report->addTitle("VarDiscover");
    o.report->addFile("VarDiscover_summary.stats");
    o.report->addFile("VarDiscover_sequins.csv");
    o.report->addFile("VarDiscover_queries.csv");
    o.report->addFile("VarDiscover_ROC.R");
    o.report->addFile("VarDiscover_LOD.R");
}