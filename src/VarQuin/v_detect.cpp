#include "tools/tools.hpp"
#include "VarQuin/v_detect.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_vcf2.hpp"

using namespace Anaquin;

typedef SequinVariant::Context Context;

extern Scripts PlotVLODR();
extern Scripts PlotVGROC();
extern Scripts PlotVCROC();
extern Scripts PlotAllele();

extern Path __output__;
extern std::string __full_command__;

inline std::string gt2str(Genotype x)
{
    switch (x)
    {
        case Genotype::Somatic:     { return "Somatic";     }
        case Genotype::Homozygous:  { return "Homozygous";  }
        case Genotype::Heterzygous: { return "Heterzygous"; }
    }
}

inline std::string ctx2Str(Context x)
{
    switch (x)
    {
        case Context::Cancer:        { return "Cancer";                    }
        case Context::LowGC:         { return "LowGC";                     }
        case Context::HighGC:        { return "HighGC";                    }
        case Context::Common:        { return "Common";                    }
        case Context::VeryLowGC:     { return "VeryLowGC";                 }
        case Context::VeryHighGC:    { return "VeryHighGC";                }
        case Context::LongHompo:     { return "LongHomopolymer";           }
        case Context::ShortHompo:    { return "ShortHomopolymer";          }
        case Context::ShortDinRep:   { return "ShortDinucleotideRepeat";   }
        case Context::LongDinRep:    { return "LongDinucleotideRepeat";    }
        case Context::ShortQuadRep:  { return "ShortQuadNucleotideRepeat"; }
        case Context::LongQuadRep:   { return "LongQuadNucleotideRepeat";  }
        case Context::ShortTrinRep:  { return "ShortTrinucleotideRepeat";  }
        case Context::LongTrinRep:   { return "LongTrinucleotideRepeat";   }
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

VDetect::EStats VDetect::analyzeE(const FileName &file, const Options &o)
{
    const auto regs = Standard::instance().r_var.regs1();
    
    VDetect::EStats stats;

    stats.v2c[Variation::SNP];
    stats.v2c[Variation::Deletion];
    stats.v2c[Variation::Inversion];
    stats.v2c[Variation::Insertion];
    stats.v2c[Variation::Duplication];
    
    if (!file.empty())
    {
        ParserVCF2::parse(file, [&](const Variant &x)
        {
            if (o.meth == VDetect::Method::Passed && x.filter != Filter::Pass)
            {
                return;
            }
            
            stats.v2c[x.type()]++;
        });
    }

    return stats;
}

VDetect::SStats VDetect::analyzeS(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VDetect::SStats stats;
    
    typedef SequinVariant::Context Context;
    
    auto muts = std::set<Variation>
    {
        Variation::SNP,
        Variation::Deletion,
        Variation::Insertion,
    };
    
    auto ctx = std::set<Context>
    {
        Context::LowGC,
        Context::HighGC,
        Context::Cancer,
        Context::Common,
        Context::LongHompo,
        Context::VeryLowGC,
        Context::VeryHighGC,
        Context::ShortDinRep,
        Context::LongDinRep,
        Context::ShortHompo,
        Context::LongQuadRep,
        Context::LongTrinRep,
        Context::ShortQuadRep,
        Context::ShortTrinRep,
    };

    for (auto &i : ctx)  { stats.g2c[i]; }
    for (auto &i : muts) { stats.v2c[i]; }
    
    o.analyze(file);

    ParserVCF2::parse(file, [&](const Variant &x)
    {
        if (o.meth == VDetect::Method::Passed && x.filter != Filter::Pass)
        {
            return;
        }
        
        const auto &r = Standard::instance().r_var;
        
        auto findMatch = [&](const Variant &query)
        {
            Match m;

            m.qry = query;
            m.var = nullptr;

            // Can we match by position?
            if ((m.var = r.findVar(query.cID, query.l)))
            {
                // Match by reference allele?
                m.ref = m.var->ref == query.ref;
                
                // Match by alternative allele?
                m.alt = m.var->alt == query.alt;
            }
            
            if (m.var)
            {
                m.rID = m.var->name;
                A_ASSERT(!m.rID.empty());
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
                        m.rID = m2->id();
                        A_ASSERT(!m.rID.empty());
                    }
                }
                catch (...) {}
            }
            
            return m;
        };
        
        const auto m = findMatch(x);
        
        // Matched if the position and alleles agree
        const auto matched = m.var && m.ref && m.alt;
        
        if (matched)
        {
            const auto key = m.var->key();
            
            stats.tps.push_back(m);
            
            const auto exp = r.findAFreq(m.var->name);
            const auto obs = m.qry.allF;
            
            A_ASSERT(!isnan(exp) && !isnan(obs));
            
            // Eg: 2821292107
            const auto id = toString(key);
            
            // Add for all variants
            stats.oa.add(id, exp, obs);
            
            // Add for mutation type
            stats.m2a[m.qry.type()].add(id, exp, obs);
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
            const auto &sv = r.findSeqVar(i.var->key());
            
            // Overall performance
            stats.oc.tp()++;
            
            // Performance for each context
            stats.g2c[sv.ctx].tp()++;
            
            // Performance by variation
            stats.v2c[i.qry.type()].tp()++;
        }
    };
    
    auto forFP = [&]()
    {
        for (auto &i : stats.fps)
        {
            // Overall performance
            stats.oc.fp()++;
            
            // Performance for each mutation
            stats.v2c[i.qry.type()].fp()++;
        }
    };

    forTP();
    forFP();

    for (auto &mut : muts)
    {
        stats.v2c[mut].nr() = r.nType(mut);
        stats.v2c[mut].nq() = stats.v2c[mut].tp() + stats.v2c[mut].fp();
        stats.v2c[mut].fn() = stats.v2c[mut].nr() - stats.v2c[mut].tp();
        stats.oc.nr() += r.nType(mut);
    }

    stats.oc.fn() = stats.oc.nr() - stats.oc.tp();
    
    for (auto &i : ctx)
    {
        stats.g2c[i].nr() = r.nContext(i);
        stats.g2c[i].nq() = stats.g2c[i].tp() + stats.g2c[i].fp();
        stats.g2c[i].fn() = stats.g2c[i].nr() - stats.g2c[i].tp();
    }
    
    A_ASSERT(stats.oc.nr() >= stats.oc.fn());
 
    for (const auto &i : r.vars())
    {
        if (!stats.findTP(i.name))
        {
            VDetect::Match m;
            
            m.var = r.findVar(i.cID, i.l);
            m.rID = i.name;
            A_ASSERT(m.var);
            
            stats.fns.push_back(m);
        }
    }
    
    return stats;
}

static void writeQuins(const FileName &file,
                       const VDetect::SStats &ss,
                       const VDetect::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Position"
                                           % "Label"
                                           % "ReadR"
                                           % "ReadV"
                                           % "Depth"
                                           % "ExpFreq"
                                           % "ObsFreq"
                                           % "Qual"
                                           % "Genotype"
                                           % "Context"
                                           % "Mutation").str());
    for (const auto &i : r.vars())
    {
        // Can we find this sequin?
        const auto isTP = ss.findTP(i.name);

        // This shouldn't fail...
        const auto &sv = r.findSeqVar(i.key());

        if (isTP)
        {
            // Called variant (if found)
            const auto &c = isTP->qry;
            
            o.writer->write((boost::format(format) % i.name
                                                   % i.cID
                                                   % i.l.start
                                                   % "TP"
                                                   % c.readR
                                                   % c.readV
                                                   % c.depth
                                                   % r.findAFreq(i.name)
                                                   % c.allF
                                                   % toString(c.qual)
                                                   % gt2str(sv.gt)
                                                   % ctx2Str(sv.ctx)
                                                   % var2str(i.type())).str());
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
                                                   % gt2str(sv.gt)
                                                   % ctx2Str(sv.ctx)
                                                   % var2str(i.type())).str());
        }
    }
    
    o.writer->close();
}

static void writeDetected(const FileName &file,
                          const VDetect::SStats &ss,
                          const VDetect::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Position"
                                           % "Label"
                                           % "ReadR"
                                           % "ReadV"
                                           % "Depth"
                                           % "ExpFreq"
                                           % "ObsFreq"
                                           % "Qual"
                                           % "Context"
                                           % "Mutation").str());

    auto f = [&](const std::vector<VDetect::Match> &x, const std::string &label)
    {
        for (const auto &i : x)
        {
            auto sID = (i.var && i.alt && i.ref ? i.var->name : "-");
            const auto ctx = sID != "-" ?  ctx2Str(r.findSeqVar(i.var->key()).ctx) : "-";

            o.writer->write((boost::format(format) % (i.rID.empty() ? "-" : i.rID)
                                                   % i.qry.cID
                                                   % i.qry.l.start
                                                   % label
                                                   % i.qry.readR
                                                   % i.qry.readV
                                                   % i.qry.depth
                                                   % (sID != "-" ? std::to_string(r.findAFreq(sID)) : "-")
                                                   % i.qry.allF
                                                   % toString(i.qry.qual)
                                                   % ctx
                                                   % var2str(i.qry.type())).str());
        }
    };
    
    f(ss.tps, "TP");
    f(ss.fps, "FP");

    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &endo,
                         const FileName &seqs,
                         const VDetect::EStats &es,
                         const VDetect::SStats &ss,
                         const VDetect::Options &o)
{
    const auto &r = Standard::instance().r_var;

    extern FileName VCFRef();
    extern FileName BedRef();
    extern FileName MixRef();

    auto germline = [&]()
    {
        const auto summary = "-------VarDetect Output Results\n\n"
                             "-------VarDetect Output\n\n"
                             "       Reference variant annotation:      %1%\n"
                             "       Reference coordinate annotation:   %2%\n\n"
                             "       User identified variants (sample): %3%\n"
                             "       User identified variants (sequin): %59%\n\n"
                             "       Number of sample variants (sequin regions): %60%\n"
                             "       Number of sequin variants (sequin regions): %61%\n\n"
                             "-------Reference variants by mutation\n\n"
                             "       SNPs:   %4%\n"
                             "       Indels: %5%\n"
                             "       Total:  %6%\n\n"
                             "-------Reference variants by context\n\n"
                             "       Common:                       %31%\n"
                             "       Very Low GC:                  %32%\n"
                             "       Low GC:                       %33%\n"
                             "       High GC:                      %34%\n"
                             "       Very High GC:                 %35%\n"
                             "       Short Dinucleotide Repeat:    %36%\n"
                             "       Long Dinucleotide Repeat:     %37%\n"
                             "       Short Homopolymer:            %38%\n"
                             "       Long Homopolymer:             %39%\n"
                             "       Short Quad Nucleotide Repeat: %40%\n"
                             "       Long Quad Nucleotide Repeat:  %41%\n"
                             "       Short Trinucleotide Repeat:   %42%\n"
                             "       Long Trinucleotide Repeat:    %43%\n\n"
                             "-------Reference variants by genotype\n\n"
                             "       Homozygosity:   %44%\n"
                             "       Heterozygosity: %45%\n\n"
                             "-------Called variants by mutation\n\n"
                             "       %7% SNPs\n"
                             "       %8% indels\n"
                             "       %9% variants\n\n"
                             "-------Diagnostic performance by mutation\n\n"
                             "       True Positive:  %10% SNPs\n"
                             "       True Positive:  %11% indels\n"
                             "       True Positive:  %12% variants\n\n"
                             "       False Positive: %13% SNPs\n"
                             "       False Positive: %14% indels\n"
                             "       False Positive: %15% variants\n\n"
                             "       False Negative: %16% SNPs\n"
                             "       False Negative: %17% indels\n"
                             "       False Negative: %18% variants\n\n"
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
                             "       FDR Rate:    %30$.4f\n\n"
                             "-------Diagnostic performance by context\n\n"
                             "       *Low GC\n"
                             "       Sensitivity: %46$.4f\n\n"
                             "       *High GC\n"
                             "       Sensitivity: %47$.4f\n\n"
                             "       *Common\n"
                             "       Sensitivity: %48$.4f\n\n"
                             "       *Long Homopolymer\n"
                             "       Sensitivity: %49$.4f\n\n"
                             "       *Very Low GC\n"
                             "       Sensitivity: %50$.4f\n\n"
                             "       *Very High GC\n"
                             "       Sensitivity: %51$.4f\n\n"
                             "       *Short Dinucleotide Repeat\n"
                             "       Sensitivity: %52$.4f\n\n"
                             "       *Long Dinucleotide Repeat\n"
                             "       Sensitivity: %53$.4f\n\n"
                             "       *Short Homopolymer\n"
                             "       Sensitivity: %54$.4f\n\n"
                             "       *Long Quad Nucleotide Repeat\n"
                             "       Sensitivity: %55$.4f\n\n"
                             "       *Long Trinucleotide Repeat\n"
                             "       Sensitivity: %56$.4f\n\n"
                             "       *Short Quad Nucleotide Repeat\n"
                             "       Sensitivity: %57$.4f\n\n"
                             "       *Short Trinucleotide Repeat\n"
                             "       Sensitivity: %58$.4f";

        #define D(x) (isnan(x) ? "-" : std::to_string(x))
        
        const auto &m2c = ss.v2c;
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

        #define CSN(x) D(ss.g2c.at(x).sn())

        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(summary) % VCFRef()                       // 1
                                                % BedRef()                       // 2
                                                % seqs                           // 3
                                                % r.nType(Variation::SNP)        // 4
                                                % (r.nType(Variation::Insertion) +
                                                   r.nType(Variation::Deletion)) // 5
                                                % (r.nType(Variation::SNP) +
                                                   r.nType(Variation::Insertion) +
                                                   r.nType(Variation::Deletion)) // 6
                                                % c_nSNP                         // 7
                                                % (c_nDel + c_nIns)              // 8
                                                % (c_nSNP + c_nDel + c_nIns)     // 9
                                                % tp_SNP                         // 10
                                                % (tp_Del + tp_Ins)              // 11
                                                % (tp_SNP + tp_Del + tp_Ins)     // 12
                                                % fp_SNP                         // 13
                                                % (fp_Del + fp_Ins)              // 14
                                                % (fp_SNP + fp_Del + fp_Ins)     // 15
                                                % fn_SNP                         // 16
                                                % (fn_Del + fn_Ins)              // 17
                                                % (fn_SNP + fn_Del + fn_Ins)     // 18
                                                % D(ss.oc.sn())                  // 19
                                                % D(ss.oc.pc())                  // 20
                                                % D(ss.oc.F1())                  // 21
                                                % D(1-ss.oc.pc())                // 22
                                                % D(snp.sn())                    // 23
                                                % D(snp.pc())                    // 24
                                                % D(snp.F1())                    // 25
                                                % D(1-snp.pc())                  // 26
                                                % D(ind.sn())                    // 27
                                                % D(ind.pc())                    // 28
                                                % D(ind.F1())                    // 29
                                                % D(1-ind.pc())                  // 30
                                                % D(r.nContext(Context::Common))       // 31
                                                % D(r.nContext(Context::VeryLowGC))    // 32
                                                % D(r.nContext(Context::LowGC))        // 33
                                                % D(r.nContext(Context::HighGC))       // 34
                                                % D(r.nContext(Context::VeryHighGC))   // 35
                                                % D(r.nContext(Context::ShortDinRep))  // 36
                                                % D(r.nContext(Context::LongDinRep))   // 37
                                                % D(r.nContext(Context::ShortHompo))   // 38
                                                % D(r.nContext(Context::LongHompo))    // 39
                                                % D(r.nContext(Context::ShortQuadRep)) // 40
                                                % D(r.nContext(Context::LongQuadRep))  // 41
                                                % D(r.nContext(Context::ShortTrinRep)) // 42
                                                % D(r.nContext(Context::LongTrinRep))  // 43
                                                % D(r.nGeno(Genotype::Homozygous))     // 44
                                                % D(r.nGeno(Genotype::Heterzygous))    // 45
                                                % CSN(Context::LowGC)                  // 46
                                                % CSN(Context::HighGC)                 // 47
                                                % CSN(Context::Common)                 // 48
                                                % CSN(Context::LongHompo)              // 49
                                                % CSN(Context::VeryLowGC)              // 50
                                                % CSN(Context::VeryHighGC)             // 51
                                                % CSN(Context::ShortDinRep)            // 52
                                                % CSN(Context::LongDinRep)             // 53
                                                % CSN(Context::ShortHompo)             // 54
                                                % CSN(Context::LongQuadRep)            // 55
                                                % CSN(Context::LongTrinRep)            // 56
                                                % CSN(Context::ShortQuadRep)           // 57
                                                % CSN(Context::ShortTrinRep)           // 58
                                                % (endo.empty() ? "-" : endo)          // 59
                                                % "????" //(endo.empty() ? "-" : toString(es.found)) // 60
                                                % (c_nSNP + c_nDel + c_nIns)           // 61
                         ).str());
    };
    
    germline();
    o.writer->close();
}

template <typename T, typename O> void writeVCF(const FileName &file, const T &x, const O &o)
{
    const auto head = "##fileformat=VCFv4.1\n"
                      "##reference=https://www.sequin.xyz\n"
                      "##INFO=<ID=AF,Number=A,Type=Float,Description=""Allele Frequency"">\n"
                      "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write(head);
    
    for (const auto &i : x)
    {
        const auto var = i.var ? i.var : &i.qry;
        const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t.\t.\tAF=%6%";
        
        o.writer->write((boost::format(format) % var->cID
                                               % var->l.start
                                               % var->name
                                               % var->ref
                                               % var->alt
                                               % var->allF).str());
    }
    
    o.writer->close();
}

void VDetect::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto es = analyzeE(endo, o);
    const auto ss = analyzeS(seqs, o);

    o.info("TP: " + std::to_string(ss.oc.tp()));
    o.info("FP: " + std::to_string(ss.oc.fp()));
    o.info("FN: " + std::to_string(ss.oc.fn()));

    o.info("Generating statistics");

    /*
     * Generating VarDetect_sequins.csv
     */
    
    writeQuins("VarDetect_sequins.csv", ss, o);

    /*
     * Generating VarDetect_summary.stats
     */
    
    writeSummary("VarDetect_summary.stats", endo, seqs, es, ss, o);
    
    /*
     * Generating VarDetect_detected.csv
     */
    
    writeDetected("VarDetect_detected.csv", ss, o);
    
    /*
     * Generating VarDetect_ROC.R
     */
    
    o.generate("VarDetect_ROC.R");
    o.writer->open("VarDetect_ROC.R");
    o.writer->write(createVGROC("VarDetect_detected.csv", "data$Depth", "'FP'"));
    o.writer->close();
    
    /*
     * Generating VarData_TP.vcf
     */
    
    writeVCF("VarDetect_TP.vcf", ss.tps, o);
    
    /*
     * Generating VarData_FP.vcf
     */
    
    writeVCF("VarData_FP.vcf", ss.fps, o);
    
    /*
     * Generating VarData_FN.vcf
     */
    
    writeVCF("VarData_FN.vcf", ss.fns, o);
}
