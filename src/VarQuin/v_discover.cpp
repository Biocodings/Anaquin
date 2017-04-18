#include "tools/tools.hpp"
#include "VarQuin/v_discover.hpp"

using namespace Anaquin;

typedef SeqVariant::Group Group;

extern Scripts PlotVLODR();
extern Scripts PlotVGROC();
extern Scripts PlotVCROC();
extern Scripts PlotAllele();

extern Path __output__;
extern std::string __full_command__;

inline std::string type2str(Variation type)
{
    switch (type)
    {
        case Variation::SNP:       { return "SNP";       }
        case Variation::Deletion:  { return "Deletion";  }
        case Variation::Insertion: { return "Insertion"; }
    }
}

inline std::string grp2Str(Group grp)
{
    switch (grp)
    {
        case Group::Cosmic:        { return "Cosmic";                    }
        case Group::LowGC:         { return "LowGC";                     }
        case Group::HighGC:        { return "HighGC";                    }
        case Group::NA12878:       { return "NA12878";                   }
        case Group::VeryLowGC:     { return "VeryLowGC";                 }
        case Group::VeryHighGC:    { return "VeryHighGC";                }
        case Group::LongHompo:     { return "LongHomopolymer";           }
        case Group::ShortHompo:    { return "ShortHomopolymer";          }
        case Group::ShortDinRep:   { return "ShortDinucleotideRepeat";   }
        case Group::LongDinRep:    { return "LongDinucleotideRepeat";    }
        case Group::ShortQuadRep:  { return "ShortQuadNucleotideRepeat"; }
        case Group::LongQuadRep:   { return "LongQuadNucleotideRepeat";  }
        case Group::ShortTrinRep:  { return "ShortTrinucleotideRepeat";  }
        case Group::LongTrinRep:   { return "LongTrinucleotideRepeat";   }
    }
}

static Scripts createVGROC(const FileName &file, const std::string &score, const std::string &refRat)
{
    return (boost::format(PlotVGROC()) % date()
                                       % __full_command__
                                       % __output__
                                       % file
                                       % score
                                       % refRat).str();
}

static Scripts createVCROC(const FileName &file, const std::string &score, const std::string &refRat)
{
    return (boost::format(PlotVCROC()) % date()
                                       % __full_command__
                                       % __output__
                                       % file
                                       % score
                                       % refRat).str();
}

VDiscover::Stats VDiscover::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VDiscover::Stats stats;
    
    typedef SeqVariant::Group Group;
    
    auto muts = std::set<Variation>
    {
        Variation::SNP,
        Variation::Deletion,
        Variation::Insertion,
    };
    
    auto grps = std::set<Group>
    {
        Group::LowGC,
        Group::HighGC,
        Group::Cosmic,
        Group::NA12878,
        Group::LongHompo,
        Group::VeryLowGC,
        Group::VeryHighGC,
        Group::ShortDinRep,
        Group::LongDinRep,
        Group::ShortHompo,
        Group::LongQuadRep,
        Group::LongTrinRep,
        Group::ShortQuadRep,
        Group::ShortTrinRep,
    };

    for (auto &i : grps) { stats.g2c[i]; }
    for (auto &i : muts) { stats.m2c[i]; }
    
    o.analyze(file);

    readVFile(file, [&](const ParserVCF::Data &x, const ParserProgress &p)
    {
        if (p.i && !(p.i % 100000))
        {
            o.wait(std::to_string(p.i));
        }
        
        const auto &r = Standard::instance().r_var;
        
        auto findMatch = [&](const ParserVCF::Data &query)
        {
            VariantMatch m;

            m.query = query;
            m.seqByPos = nullptr;
            
            // Can we match by position?
            if ((m.seqByPos = r.findVar(query.cID, query.l)))
            {
                // Match by reference allele?
                m.ref = m.seqByPos->ref == query.ref;
                
                // Match by alternative allele?
                m.alt = m.seqByPos->alt == query.alt;
            }
            
            if (m.seqByPos)
            {
                m.rReg = m.seqByPos->name;
                A_ASSERT(!m.rReg.empty());
            }
            else
            {
                MergedIntervals<> inters;
                
                try
                {
                    // We should to search where the FPs are
                    inters = r.mInters(x.cID);
                    
                    A_ASSERT(inters.size());
                    
                    const auto m2 = inters.contains(x.l);
                    
                    // Can we find the corresponding region for the FP?
                    if (m2)
                    {
                        m.rReg = m2->id();
                        A_ASSERT(!m.rReg.empty());
                    }
                }
                catch (...) {}
            }
            
            return m;
        };
        
        const auto m = findMatch(x);
        
        // Matched if the position and alleles agree
        const auto matched = m.seqByPos && m.ref && m.alt;
        
        if (matched)
        {
            const auto key = m.seqByPos->key();
            
            stats.tps.push_back(m);
            
            const auto exp = r.findAFreq(m.seqByPos->name);
            const auto obs = m.query.allF;
            
            A_ASSERT(!isnan(exp) && !isnan(obs));
            
            // Eg: 2821292107
            const auto id = toString(key);
            
            // Add for all variants
            stats.oa.add(id, exp, obs);
            
            // Add for mutation type
            stats.m2a[m.query.type()].add(id, exp, obs);
        }
        else
        {
            // FP because the variant is not found in the reference
            stats.fps.push_back(m);
        }
    });
    
    o.info("Aggregating statistics");

    /*
     * Determine the quantification limits
     */
    
    stats.oa.limit = stats.oa.limitQuant();

    for (auto &i : stats.m2a)
    {
        i.second.limit = i.second.limitQuant();
    }

    /*
     * Determine the classification performance
     */

    auto forTP = [&]()
    {
        for (auto &i : stats.tps)
        {
            // This shouldn't fail...
            const auto &sv = r.findSeqVar(i.seqByPos->key());
            
            // Overall performance
            stats.oc.tp()++;
            
            // Performance for each group
            stats.g2c[sv.group].tp()++;
            
            // Performance for each mutation
            stats.m2c[i.query.type()].tp()++;
        }
    };
    
    auto forFP = [&]()
    {
        for (auto &i : stats.fps)
        {
            // Overall performance
            stats.oc.fp()++;
            
            // Performance for each mutation
            stats.m2c[i.query.type()].fp()++;
        }
    };

    forTP();
    forFP();

    for (auto &mut : muts)
    {
        stats.m2c[mut].nr() = r.nType(mut);
        stats.m2c[mut].nq() = stats.m2c[mut].tp() + stats.m2c[mut].fp();
        stats.m2c[mut].fn() = stats.m2c[mut].nr() - stats.m2c[mut].tp();
        stats.oc.nr() += r.nType(mut);
    }

    stats.oc.fn() = stats.oc.nr() - stats.oc.tp();
    
    for (auto &grp : grps)
    {
        stats.g2c[grp].nr() = r.nGroup(grp);
        stats.g2c[grp].nq() = stats.g2c[grp].tp() + stats.g2c[grp].fp();
        stats.g2c[grp].fn() = stats.g2c[grp].nr() - stats.g2c[grp].tp();
    }
    
    A_ASSERT(stats.oc.nr() >= stats.oc.fn());
    
    return stats;
}

static void writeQuins(const FileName &file,
                       const VDiscover::Stats &stats,
                       const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "ChrID"
                                           % "Position"
                                           % "Label"
                                           % "ReadR"
                                           % "ReadV"
                                           % "Depth"
                                           % "ExpFreq"
                                           % "ObsFreq"
                                           % "Pval"
                                           % "Qual"
                                           % "Group"
                                           % "Type").str());
    for (const auto &i : r.vars())
    {
        // Can we find this sequin?
        const auto isTP = stats.findTP(i.name);

        // This shouldn't fail...
        const auto &sv = r.findSeqVar(i.key());

        if (isTP)
        {
            // Called variant (if found)
            const auto &c = isTP->query;
            
            o.writer->write((boost::format(format) % i.name
                                                   % i.cID
                                                   % i.l.start
                                                   % "TP"
                                                   % c.readR
                                                   % c.readV
                                                   % c.depth
                                                   % r.findAFreq(i.name)
                                                   % c.allF
                                                   % ld2ss(c.p)
                                                   % toString(c.qual)
                                                   % grp2Str(sv.group)
                                                   % type2str(i.type())).str());
        }

        // Failed to detect the variant
        else
        {
            o.writer->write((boost::format(format) % i.name
                                                   % i.cID
                                                   % i.l.start
                                                   % "FN"
                                                   % "-"
                                                   % "-"
                                                   % "-"
                                                   % r.findAFreq(i.name)
                                                   % "-"
                                                   % "-"
                                                   % "-"
                                                   % grp2Str(sv.group)
                                                   % type2str(i.type())).str());
        }
    }
    
    o.writer->close();
}

static void writeDetected(const FileName &file, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "ChrID"
                                           % "Position"
                                           % "Label"
                                           % "ReadR"
                                           % "ReadV"
                                           % "Depth"
                                           % "ExpFreq"
                                           % "ObsFreq"
                                           % "Pval"
                                           % "Qual"
                                           % "Group"
                                           % "Type").str());

    auto f = [&](const std::vector<VariantMatch> &x, const std::string &label)
    {
        for (const auto &i : x)
        {
            auto sID = (i.seqByPos && i.alt && i.ref ? i.seqByPos->name : "-");
            const auto grp = sID != "-" ?  grp2Str(r.findSeqVar(i.seqByPos->key()).group) : "-";

            o.writer->write((boost::format(format) % (i.rReg.empty() ? "-" : i.rReg)
                                                   % i.query.cID
                                                   % i.query.l.start
                                                   % label
                                                   % i.query.readR
                                                   % i.query.readV
                                                   % i.query.depth
                                                   % (sID != "-" ? std::to_string(r.findAFreq(sID)) : "-")
                                                   % i.query.allF
                                                   % ld2ss(i.query.p)
                                                   % toString(i.query.qual)
                                                   % grp
                                                   % type2str(i.query.type())).str());
        }
    };
    
    f(stats.tps, "TP");
    f(stats.fps, "FP");

    o.writer->close();
}

static void writeSummary(const FileName &file, const FileName &src, const VDiscover::Stats &stats, const VDiscover::Options &o)
{
    const auto &r = Standard::instance().r_var;

    extern FileName VCFRef();
    extern FileName BedRef();
    extern FileName MixRef();

    const auto &ss = stats;

    auto germline = [&]()
    {
        const auto summary = "-------VarDiscover Output Results\n\n"
                             "-------VarDiscover Output\n\n"
                             "       Reference variant annotations:    %1%\n"
                             "       Reference coordinate annotations: %2%\n"
                             "       User identified variants:         %3%\n\n"
                             "-------Reference variants by type\n\n"
                             "       %4% SNPs\n"
                             "       %5% indels\n"
                             "       %6% variants\n\n"
                             "-------Called variants by type\n\n"
                             "       %7% SNPs\n"
                             "       %8% indels\n"
                             "       %9% variants\n\n"
                             "-------Identification of synthetic variants\n\n"
                             "       True Positive:  %10% SNPS\n"
                             "       True Positive:  %11% indels\n"
                             "       True Positive:  %12% variants\n\n"
                             "       False Positive: %13% SNPs\n"
                             "       False Positive: %14% indels\n"
                             "       False Positive: %15% variants\n\n"
                             "       False Negative: %16% SNPs\n"
                             "       False Negative: %17% indels\n"
                             "       False Negative: %18% variants\n\n"
                             "-------Diagnostic Performance (Synthetic)\n\n"
                             "       *Variants\n"
                             "       Sensitivity: %19$.4f\n"
                             "       Precision:   %20$.4f\n"
                             "       F1 Score:    %21$.4f\n"
                             "       FDR Rate:    %22$.4f\n\n"
                             "       *SNPs\n"
                             "       Sensitivity: %23$.4f\n"
                             "       Precision:   %24$.4f\n"
                             "       F1 Score:    %25$.4f\n"
                             "       FDR Rate:    %26$.4f\n\n"
                             "       *Indels\n"
                             "       Sensitivity: %27$.4f\n"
                             "       Precision:   %28$.4f\n"
                             "       F1 Score:    %29$.4f\n"
                             "       FDR Rate:    %30$.4f";

        #define D(x) (isnan(x) ? "-" : std::to_string(x))
        
        const auto &m2c = stats.m2c;
        const auto &snp = m2c.at(Variation::SNP);
        const auto &del = m2c.at(Variation::Deletion);
        const auto &ins = m2c.at(Variation::Insertion);
        
        const auto c_nSNP = snp.nq();
        const auto c_nDel = del.nq();
        const auto c_nIns = ins.nq();
        
        const auto tp_SNP = snp.tp();
        const auto tp_Del = del.tp();
        const auto tp_Ins = ins.tp();

        const auto fp_SNP = snp.fp();
        const auto fp_Del = del.fp();
        const auto fp_Ins = ins.fp();
        
        const auto fn_SNP = snp.fn();
        const auto fn_Del = del.fn();
        const auto fn_Ins = ins.fn();
        
        auto ind = del;
        ind += ins;

        o.generate(file);
        o.writer->open("VarDiscover_summary.stats");
        o.writer->write((boost::format(summary) % VCFRef()                      // 1
                                                % BedRef()                      // 2
                                                % src                           // 3
                                                % r.countSNP()                  // 4
                                                % r.countInd()                  // 5
                                                % (r.countSNP() + r.countInd()) // 6
                                                % c_nSNP                        // 7
                                                % (c_nDel + c_nIns)             // 8
                                                % (c_nSNP + c_nDel + c_nIns)    // 9
                                                % tp_SNP                        // 10
                                                % (tp_Del + tp_Ins)             // 11
                                                % (tp_SNP + tp_Del + tp_Ins)    // 12
                                                % fp_SNP                        // 13
                                                % (fp_Del + fp_Ins)             // 14
                                                % (fp_SNP + fp_Del + fp_Ins)    // 15
                                                % fn_SNP                        // 16
                                                % (fn_Del + fn_Ins)             // 17
                                                % (fn_SNP + fn_Del + fn_Ins)    // 18
                                                % D(stats.oc.sn())              // 19
                                                % D(stats.oc.pc())              // 20
                                                % D(stats.oc.F1())              // 21
                                                % D(1-stats.oc.pc())            // 22
                                                % D(snp.sn())                   // 23
                                                % D(snp.pc())                   // 24
                                                % D(snp.F1())                   // 25
                                                % D(1-snp.pc())                 // 26
                                                % D(ind.sn())                   // 27
                                                % D(ind.pc())                   // 28
                                                % D(ind.F1())                   // 29
                                                % D(1-ind.pc())                 // 30
                         ).str());
    };
    
    auto somatic = [&]()
    {
        const auto lm = ss.oa.linear(true);

        const auto summary = "-------VarDiscover Output Results\n\n"
                             "-------VarDiscover Output\n\n"
                             "       Reference variant annotations:    %1%\n"
                             "       Reference coordinate annotations: %2%\n"
                             "       User identified variants:         %3%\n\n"
                             "-------Reference variant annotations\n\n"
                             "       Synthetic: %4% SNPs\n"
                             "       Synthetic: %5% indels\n"
                             "       Synthetic: %6% variants\n\n"
                             "-------User identified variants\n\n"
                             "       Synthetic: %7% SNPs\n"
                             "       Synthetic: %8% indels\n"
                             "       Synthetic: %9% variants\n\n"
                             "       Detection Sensitivity: %10% (attomol/ul) (%11%)\n\n"
                             "-------Identification of synthetic variants\n\n"
                             "       True Positive:  %12% SNPS\n"
                             "       True Positive:  %13% indels\n"
                             "       True Positive:  %14% variants\n\n"
                             "       False Positive: %15% SNPs\n"
                             "       False Positive: %16% indels\n"
                             "       False Positive: %17% variants\n\n"
                             "       False Negative: %18% SNPs\n"
                             "       False Negative: %19% indels\n"
                             "       False Negative: %20% variants\n\n"
                             "-------Diagnostic Performance (Synthetic)\n\n"
                             "       *Variants\n"
                             "       Sensitivity: %21$.4f\n"
                             "       Precision:   %22$.4f\n"
                             "       FDR Rate:    %23$.4f\n\n"
                             "       *SNPs\n"
                             "       Sensitivity: %24$.4f\n"
                             "       Precision:   %25$.4f\n"
                             "       FDR Rate:    %26$.4f\n\n"
                             "       *Indels\n"
                             "       Sensitivity: %27$.4f\n"
                             "       Precision:   %28$.4f\n"
                             "       FDR Rate:    %29$.4f\n\n"
                             "-------Overall linear regression (log2 scale)\n\n"
                             "      Slope:       %30%\n"
                             "      Correlation: %31%\n"
                             "      R2:          %32%\n"
                             "      F-statistic: %33%\n"
                             "      P-value:     %34%\n"
                             "      SSM:         %35%, DF: %36%\n"
                             "      SSE:         %37%, DF: %38%\n"
                             "      SST:         %39%, DF: %40%\n";
        o.generate(file);
        o.writer->open("VarDiscover_summary.stats");
        o.writer->write((boost::format(summary) % VCFRef()                   // 1
                                                % BedRef()                   // 2
                                                % src                        // 3
                                                % r.countSNP()            // 4
                                                % r.countInd()            // 5
                                                % (r.countSNP() + r.countInd()) // 6
                                                % "??" //r.countSNPGen()            // 7
                                                % "??" //r.countIndGen()            // 8
                                                % "??" //(r.countSNPGen() + r.countIndGen()) // 9
                                                % "??" //ss.vData.countSNPSyn()  // 10
                                                % "??" //ss.vData.countIndSyn()  // 11
                                                % "??" //ss.vData.countVarSyn()  // 12
                                                % "??" //ss.vars.limit.abund     // 13
                                                % "??" //ss.vars.limit.id        // 14
                                                % "??" //ss.countSNP_TP_Syn()    // 15
                                                % "??" //ss.countInd_TP_Syn()    // 16
                                                % "??" //ss.countVar_TP_Syn()    // 17
                                                % "??" //ss.countSNP_FP_Syn()    // 18
                                                % "??" //ss.countInd_FP_Syn()    // 19
                                                % "??" //ss.countVar_FP_Syn()    // 20
                                                % "??" //ss.countSNP_FnSyn()     // 21
                                                % "??" //ss.countInd_FnSyn()     // 22
                                                % "??" //ss.countVar_FnSyn()     // 23
                                                % "??" //ss.countVarSnSyn()      // 24
                                                % "??" //ss.countVarPC_Syn()     // 25
                                                % "??" //(1-ss.countVarPC_Syn()) // 26
                                                % "??" //ss.countSNPSnSyn()      // 27
                                                % "??" //ss.countSNPPC_Syn()     // 28
                                                % "??" //(1-ss.countSNPPC_Syn()) // 29
                                                % lm.m                       // 30
                                                % lm.r                       // 31
                                                % lm.R2                      // 32
                                                % lm.F                       // 33
                                                % lm.p                       // 34
                                                % lm.SSM                     // 35
                                                % lm.SSM_D                   // 36
                                                % lm.SSE                     // 37
                                                % lm.SSE_D                   // 38
                                                % lm.SST                     // 39
                                                % lm.SST_D                   // 40
                         ).str());
    };
    
    if (r.isGermlineRef())
    {
        germline();
    }
    else
    {
        somatic();
    }

    o.writer->close();
}

void VDiscover::report(const FileName &seqs, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto ss = analyze(seqs, o);
    
    o.info("TP: " + std::to_string(ss.oc.tp()));
    o.info("FP: " + std::to_string(ss.oc.fp()));
    o.info("FN: " + std::to_string(ss.oc.fn()));

    o.info("Generating statistics");

    /*
     * Generating VarDiscover_sequins.csv
     */
    
    writeQuins("VarDiscover_sequins.csv", ss, o);

    /*
     * Generating VarDiscover_summary.stats
     */
    
    writeSummary("VarDiscover_summary.stats", seqs, ss, o);
    
    /*
     * Generating VarDiscover_detected.csv
     */
    
    writeDetected("VarDiscover_detected.csv", ss, o);
    
    /*
     * Generating VarDiscover_ROC.R
     */
    
    o.generate("VarDiscover_ROC.R");
    o.writer->open("VarDiscover_ROC.R");
    
//    if (__countP__ >= __countD__)
//    {
//        o.info("P-value for scoring");
//        o.writer->write(createVCROC("VarDiscover_detected.csv", "1-data$Pval", "-1"));
//    }
//    else
//    {
        o.info("Depth for scoring");
        o.writer->write(createVGROC("VarDiscover_detected.csv", "data$Depth", "'FP'"));
//    }
    
    o.writer->close();
    
    if (!r.isGermlineRef())
    {
        /*
         * Generating VarDiscover_allele.R
         */
        
        o.generate("VarDiscover_allele.R");
        o.writer->open("VarDiscover_allele.R");
        o.writer->write(RWriter::createRLinear("VarDiscover_sequins.csv",
                                               o.work,
                                               "Allele Frequency",
                                               "Expected Allele Frequency (log2)",
                                               "Measured Allele Frequency (log2)",
                                               "log2(data$ExpFreq)",
                                               "log2(data$ObsFreq)",
                                               "input",
                                                true,
                                                PlotAllele()));
        o.writer->close();

        /*
         * Generating VarDiscover_LODR.R
         */
        
        o.generate("VarDiscover_LODR.R");
        o.writer->open("VarDiscover_LODR.R");
        o.writer->write(RWriter::createScript("VarDiscover_detected.csv", PlotVLODR()));
        o.writer->close();
    }
}
