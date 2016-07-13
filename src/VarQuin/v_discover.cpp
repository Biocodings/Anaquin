#include "data/convert.hpp"
#include "VarQuin/v_discover.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVLOD();

// Defined in resources.cpp
extern Scripts PlotVROC();

// Defined in resources.cpp
extern Scripts PlotGermlineROC();

static Counts __countD__ = 0;
static Counts __countP__ = 0;

struct VDiscoverImpl : public VCFDataUser
{
    VDiscover::Stats *stats;
    const VDiscover::Options *o;

    void variantProcessed(const ParserVCF::Data &x, const ParserProgress &p) override
    {
        if (p.i && !(p.i % 100000))
        {
            o->wait(std::to_string(p.i));
        }
        
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
        
        // Always work on the queries
        stats->query[cID].af.insert(m.query.alleleFreq());

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

    o.analyze(file);

    VDiscoverImpl impl;
    impl.o = &o;
    impl.stats = &stats;

    // Read the input variant file and process them
    stats.vData = vcfData(file, o.format, &impl);

    o.info("Aggregating statistics");

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
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Pos"
                                           % "Label"
                                           % "Pval"
                                           % "ReadR"
                                           % "ReadV"
                                           % "Depth"
                                           % "ERef"
                                           % "EVar"
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
                
                // Eg: "SNP"
                const auto type = type2str(m->type());

                // Unique ID for the variant
                const auto id = (m->id + "_" + std::to_string(m->l.start) + "_" + type);

                /*
                 * Now we need to know the label for this reference variant
                 */
                
                auto f = [&](const std::map<long, VariantMatch> &x, const std::string &label)
                {
                    if (x.count(key))
                    {
                        const auto &t = x.at(key);
                        
                        o.writer->write((boost::format(format) % id
                                                               % m->l.start
                                                               % label
                                                               % (isnan(t.query.p) ? "-" : ld2ss(t.query.p))
                                                               % t.query.readR
                                                               % t.query.readV
                                                               % t.query.depth
                                                               % r.findRCon(m->id)
                                                               % r.findVCon(m->id)
                                                               % r.findAFreq(m->id)
                                                               % t.query.alleleFreq()
                                                               % type).str());
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
                
                // Eg: "SNP"
                const auto type = type2str(m->type());
                
                // Unique ID for the variant
                const auto id = (m->id + "_" + std::to_string(m->l.start) + "_" + type);
                
                o.writer->write((boost::format(format) % id
                                                       % m->l.start
                                                       % "FN"
                                                       % "NA"
                                                       % "NA"
                                                       % "NA"
                                                       % "NA"
                                                       % r.findRCon(m->id)
                                                       % r.findVCon(m->id)
                                                       % r.findAFreq(m->id)
                                                       % "NA"
                                                       % type).str());
            }
        }
    }
    
    o.writer->close();
}

static void writeQueries(const FileName &file, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Pos"
                                           % "Label"
                                           % "Ref"
                                           % "Var"
                                           % "Depth"
                                           % "ERef"
                                           % "EVar"
                                           % "EFreq"
                                           % "MFreq"
                                           % "Pval"
                                           % "Qual"
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
                    
                    if (label == "FP")
                    {
                        const auto m = inters.contains(i.query.l);
                        
                        // Can we find the corresponding region for the FP?
                        if (m)
                        {
                            sID = m->id();

                            // It has to be sequin name (eg: D_3_12)
                            assert(!sID.empty());
                        }
                    }
                    
                    const auto eRef  = sID != "-" ? r.findRCon(sID)  : NAN;
                    const auto eVar  = sID != "-" ? r.findVCon(sID)  : NAN;
                    const auto eFreq = sID != "-" ? r.findAFreq(sID) : NAN;

                    const auto pval  = (isnan(i.query.p) ? "-" : ld2ss(i.query.p));
                    const auto mFreq = i.query.alleleFreq();
                    
                    o.writer->write((boost::format(format) % sID
                                                           % i.query.l.start
                                                           % label
                                                           % i.query.readR
                                                           % i.query.readV
                                                           % i.query.depth
                                                           % eRef
                                                           % eVar
                                                           % eFreq
                                                           % mFreq
                                                           % pval
                                                           % i.query.qual
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

    const auto lm = stats.vars.linear(true);
    
    // Calcluate the quantification point with logarithm
    auto ms = stats.vars.limitQuant(true);
    
    // Remember the break-point is on the log2-scale, we'll need to convert it back
    ms.b = pow(2, ms.b);
    
    Counts n_below = 0;
    Counts n_above = 0;
    
    // For each query chromosome...
    for (const auto &i : stats.query)
    {
        // For each genomic chromosome...
        if (!Standard::isSynthetic(i.first))
        {
            for (const auto &j : i.second.af)
            {
                if (j >= ms.b)
                {
                    n_above++;
                }
                else
                {
                    n_below++;
                }
            }
        }
    }

    // For each reference chromosome...
    for (const auto &i : stats.data)
    {
        // For each genomic chromosome...
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
    
    const auto mm = r.findVar(ChrIS, stol(ms.id));
    assert(mm);
    
    auto writeNoScale = [&]()
    {
        const auto summary = "-------VarDiscover Output Results\n\n"
                             "-------VarDiscover Output\n\n"
                             "       Reference variant annotations:    %1%\n"
                             "       Reference coordinate annotations: %2%\n"
                             "       Sequin mixture file:              %3%\n"
                             "       User identified variants:         %4%\n\n"
                             "-------Reference annotated variants\n\n"
                             "       Synthetic: %5% SNPs\n"
                             "       Synthetic: %6% indels\n"
                             "       Synthetic: %7% variants\n\n"
                             "-------User identified variants\n\n"
                             "       Synthetic: %8% SNPs\n"
                             "       Synthetic: %9% indels\n"
                             "       Synthetic: %10% variants\n\n"
                             "-------Identification of synthetic variants\n\n"
                             "       True Positive:  %11% SNPS\n"
                             "       True Positive:  %12% indels\n"
                             "       True Positive:  %13% variants\n\n"
                             "       False Positive: %14% SNPs\n"
                             "       False Positive: %15% indels\n"
                             "       False Positive: %16% variants\n\n"
                             "       False Negative: %17% SNPs\n"
                             "       False Negative: %19% indels\n"
                             "       False Negative: %19% variants\n\n"
                             "-------Diagnostic Performance (Synthetic)\n\n"
                             "       *Variants\n"
                             "       Sensitivity: %20$.4f\n"
                             "       Precision:   %21$.4f\n"
                             "       FDR Rate:    %22$.4f\n\n"
                             "       *SNVs\n"
                             "       Sensitivity: %23$.4f\n"
                             "       Precision:   %24$.4f\n"
                             "       FDR Rate:    %25$.4f\n\n"
                             "       *Indels\n"
                             "       Sensitivity: %26$.4f\n"
                             "       Precision:   %27$.4f\n"
                             "       FDR Rate:    %28$.4f\n";
        o.generate(file);
        o.writer->open("VarDiscover_summary.stats");
        o.writer->write((boost::format(summary) % VCFRef()                   // 1
                                                % BedRef()                   // 2
                                                % MixRef()                   // 3
                                                % src                        // 4
                                                % r.countSNPSyn()            // 5
                                                % r.countIndSyn()            // 6
                                                % (r.countSNPSyn() + r.countIndSyn())
                                                % stats.vData.countSNPSyn()  // 8
                                                % stats.vData.countIndSyn()  // 9
                                                % stats.vData.countVarSyn()  // 10
                                                % stats.countSNP_TP_Syn()    // 11
                                                % stats.countInd_TP_Syn()    // 12
                                                % stats.countVar_TP_Syn()    // 13
                                                % stats.countSNP_FP_Syn()    // 14
                                                % stats.countInd_FP_Syn()    // 15
                                                % stats.countVar_FP_Syn()    // 16
                                                % stats.countSNP_FN_Syn()    // 17
                                                % stats.countInd_FN_Syn()    // 18
                                                % stats.countVar_FN_Syn()    // 19
                                                % stats.countVarSN_Syn()     // 20
                                                % stats.countVarPC_Syn()     // 21
                                                % (1-stats.countVarPC_Syn()) // 22
                                                % stats.countSNPSN_Syn()     // 23
                                                % stats.countSNPPC_Syn()     // 24
                                                % (1-stats.countSNPPC_Syn()) // 25
                                                % stats.countIndSN_Syn()     // 26
                                                % stats.countIndPC_Syn()     // 27
                                                % (1-stats.countIndPC_Syn()) // 28
                         ).str());
    };
    
    auto writeScale = [&]()
    {
        const auto summary = "-------VarDiscover Output Results\n\n"
                             "-------VarDiscover Output\n\n"
                             "       Reference variant annotations:    %1%\n"
                             "       Reference coordinate annotations: %2%\n"
                             "       Sequin mixture file:              %3%\n"
                             "       User identified variants:         %4%\n\n"
                             "-------Reference variant annotations\n\n"
                             "       Synthetic: %5% SNPs\n"
                             "       Synthetic: %6% indels\n"
                             "       Synthetic: %7% variants\n\n"
                             "-------User identified variants\n\n"
                             "       Synthetic: %11% SNPs\n"
                             "       Synthetic: %12% indels\n"
                             "       Synthetic: %13% variants\n\n"
                             "       Detection Sensitivity: %14% (attomol/ul) (%15%)\n\n"
                             "-------Identification of synthetic variants\n\n"
                             "       True Positive:  %16% SNPS\n"
                             "       True Positive:  %17% indels\n"
                             "       True Positive:  %18% variants\n\n"
                             "       False Positive: %19% SNPs\n"
                             "       False Positive: %20% indels\n"
                             "       False Positive: %21% variants\n\n"
                             "       False Negative: %22% SNPs\n"
                             "       False Negative: %23% indels\n"
                             "       False Negative: %24% variants\n\n"
                             "-------Diagnostic Performance (Synthetic)\n\n"
                             "       *Variants\n"
                             "       Sensitivity: %25$.4f\n"
                             "       Precision:   %26$.4f\n"
                             "       FDR Rate:    %27$.4f\n\n"
                             "       *SNVs\n"
                             "       Sensitivity: %28$.4f\n"
                             "       Precision:   %29$.4f\n"
                             "       FDR Rate:    %30$.4f\n\n"
                             "       *Indels\n"
                             "       Sensitivity: %31$.4f\n"
                             "       Precision:   %32$.4f\n"
                             "       FDR Rate:    %33$.4f\n\n"
                             "-------Limit of Quantification (LOQ)\n"
                             "      *Estimated by piecewise segmented regression\n\n"
                             "       Break LOQ: %34% attomol/ul (%35%)\n\n"
                             "      *Below LOQ\n"
                             "       Intercept:   %36%\n"
                             "       Slope:       %37%\n"
                             "       Correlation: %38%\n"
                             "       R2:          %39%\n"
                             "       Genome:      %40%\n\n"
                             "      *Above LOQ\n"
                             "       Intercept:   %41%\n"
                             "       Slope:       %42%\n"
                             "       Correlation: %43%\n"
                             "       R2:          %44%\n"
                             "       Genome:      %45%\n\n"
                             "-------Overall linear regression (log2 scale)\n\n"
                             "      Slope:       %46%\n"
                             "      Correlation: %47%\n"
                             "      R2:          %48%\n"
                             "      F-statistic: %49%\n"
                             "      P-value:     %50%\n"
                             "      SSM:         %51%, DF: %52%\n"
                             "      SSE:         %53%, DF: %54%\n"
                             "      SST:         %55%, DF: %56%\n";
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
                                                % stats.countSNP_TP_Syn()    // 16
                                                % stats.countInd_TP_Syn()    // 17
                                                % stats.countVar_TP_Syn()    // 18
                                                % stats.countSNP_FP_Syn()    // 19
                                                % stats.countInd_FP_Syn()    // 20
                                                % stats.countVar_FP_Syn()    // 21
                                                % stats.countSNP_FN_Syn()    // 22
                                                % stats.countInd_FN_Syn()    // 23
                                                % stats.countVar_FN_Syn()    // 24
                                                % stats.countVarSN_Syn()     // 25
                                                % stats.countVarPC_Syn()     // 26
                                                % (1-stats.countVarPC_Syn()) // 27
                                                % stats.countSNPSN_Syn()     // 28
                                                % stats.countSNPPC_Syn()     // 29
                                                % (1-stats.countSNPPC_Syn()) // 30
                                                % stats.countIndSN_Syn()     // 31
                                                % stats.countIndPC_Syn()     // 32
                                                % (1-stats.countIndPC_Syn()) // 33
                                                % ms.b                       // 34
                                                % mm->id                     // 35
                                                % ms.lInt                    // 36
                                                % ms.lSl                     // 37
                                                % ms.lr                      // 38
                                                % ms.lR2                     // 39
                                                % n_above                    // 40
                                                % ms.rInt                    // 41
                                                % ms.rSl                     // 42
                                                % ms.rr                      // 43
                                                % ms.rR2                     // 43
                                                % n_below                    // 45
                                                % lm.m                       // 47
                                                % lm.r                       // 46
                                                % lm.R2                      // 48
                                                % lm.F                       // 49
                                                % lm.p                       // 50
                                                % lm.SSM                     // 51
                                                % lm.SSM_D                   // 52
                                                % lm.SSE                     // 53
                                                % lm.SSE_D                   // 54
                                                % lm.SST                     // 55
                                                % lm.SST_D                   // 56
                         ).str());
    };
    
    if (r.isGermline())
    {
        writeNoScale();
    }
    else
    {
        writeScale();
    }

    o.writer->close();
}

void VDiscover::report(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

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
    
    if (!r.isGermline())
    {
        /*
         * Generating VarDiscover_allele.R
         */
        
        o.generate("VarDiscover_allele.R");
        o.writer->open("VarDiscover_allele.R");
        o.writer->write(RWriter::createScatterNeedLog("VarDiscover_sequins.csv",
                                                      "Allele Frequency",
                                                      "Expected allele frequency (log2)",
                                                      "Measured allele frequency (log2)",
                                                      "EFreq",
                                                      "MFreq",
                                                      "expected",
                                                      true));
        o.writer->close();
    }
    
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
    
    if (!r.isGermline())
    {
        o.report->addFile("VarDiscover_allele.R");
    }
}