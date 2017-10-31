#include "VarQuin/v_germ.hpp"
#include "writers/vcf_writer.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

typedef SequinVariant::Context Context;

inline bool isGerm(const Variant &x)
{
    return Standard::instance().r_var.ctx2(x) != Context::Cancer;
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

static Scripts createROC(const FileName &file, const std::string &score, const std::string &refRat)
{
    extern Scripts PlotVGROC();
    extern Path __output__;
    extern std::string __full_command__;
    
    return (boost::format(PlotVGROC()) % date()
                                       % __full_command__
                                       % __output__
                                       % file
                                       % score
                                       % refRat).str();
}

VGerm::EStats VGerm::analyzeE(const FileName &file, const Options &o)
{
    const auto r2 = Standard::instance().r_var.regs2();
    
    VGerm::EStats stats;

    stats.g2c[Genotype::Homozygous];
    stats.g2c[Genotype::Heterzygous];
    stats.v2c[Variation::SNP];
    stats.v2c[Variation::Deletion];
    stats.v2c[Variation::Inversion];
    stats.v2c[Variation::Insertion];
    stats.v2c[Variation::Duplication];
    
    if (!file.empty())
    {
        ParserVCF::parse(file, [&](const Variant &x)
        {
            if (o.meth == VGerm::Method::Passed && x.filter != Filter::Pass)
            {
                return;
            }
            else if (contains(r2, x.cID, x.l))
            {
                stats.g2c[x.gt]++;
                stats.v2c[x.type()]++;
            }
        });
    }

    return stats;
}

VGerm::SStats VGerm::analyzeS(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto r2 = Standard::instance().r_var.regs2();

    VGerm::SStats stats;
    
    typedef SequinVariant::Context Context;
    
    auto gts = std::set<Genotype>
    {
        Genotype::Homozygous,
        Genotype::Heterzygous
    };

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

    for (auto &i : gts)  { stats.g2c[i]; }
    for (auto &i : ctx)  { stats.c2c[i]; }
    for (auto &i : muts) { stats.v2c[i]; }
    
    o.analyze(file);

    auto wTP = VCFWriter(); wTP.open(o.work + "/VarGermline_TP.vcf");
    auto wFP = VCFWriter(); wFP.open(o.work + "/VarGermline_FP.vcf");
    
    ParserVCF::parse(file, [&](const Variant &x)
    {
        if (o.meth == VGerm::Method::Passed && x.filter != Filter::Pass)
        {
            return;
        }
        else if (!contains(r2, x.cID, x.l))
        {
            return;
        }
        
        auto findMatch = [&](const Variant &query)
        {
            Match m;

            m.qry = query;
            m.var = nullptr;

            // Can we match by position?
            if ((m.var = r.findV1(query.cID, query.l)))
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
            if (isGerm(*m.var))
            {
                wTP.write(x.hdr, x.line);
                
                const auto key = m.var->key();
                
                stats.tps.push_back(m);
                
                const auto exp = r.af(m.var->name);
                const auto obs = m.qry.allF;
                
                A_ASSERT(!isnan(exp));
                
                // Eg: 2821292107
                const auto id = toString(key);
                
                // Add for all variants
                stats.oa.add(id, exp, obs);
                
                // Add for mutation type
                stats.m2a[m.qry.type()].add(id, exp, obs);
            }
        }
        else
        {
            // Ignore anything that is somatic
            if (!r.findV2(x.cID, x.l))
            {
                wFP.write(x.hdr, x.line);
                
                // FP because the variant is not found in the reference
                stats.fps.push_back(m);
            }
        }
    });
    
    wTP.close();
    wFP.close();
    
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
            const auto &sv = r.findSeqVar1(i.var->key());
            
            // Overall performance
            stats.oc.tp()++;
            
            // Performance by genotype
            stats.g2c[sv.gt].tp()++;
            
            // Performance by context
            stats.c2c[sv.ctx].tp()++;
            
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
            
            // Performance by mutation
            stats.v2c[i.qry.type()].fp()++;
            
            // Performance by genotype
            stats.g2c[i.qry.gt].fp()++;
        }
    };

    forTP();
    forFP();

    for (auto &mut : muts)
    {
        stats.v2c[mut].nr() = r.nType1(mut);
        stats.v2c[mut].nq() = stats.v2c[mut].tp() + stats.v2c[mut].fp();
        stats.v2c[mut].fn() = stats.v2c[mut].nr() - stats.v2c[mut].tp();
        stats.oc.nr() += r.nType1(mut);
    }

    stats.oc.fn() = stats.oc.nr() - stats.oc.tp();
    
    /*
     * Performance by context
     */
    
    for (auto &i : ctx)
    {
        stats.c2c[i].nr() = r.nCtx1(i);
        stats.c2c[i].nq() = stats.c2c[i].tp() + stats.c2c[i].fp();
        stats.c2c[i].fn() = stats.c2c[i].nr() - stats.c2c[i].tp();
    }
    
    /*
     * Performance by genotype
     */
    
    for (auto &i : gts)
    {
        stats.g2c[i].nr() = r.nGeno1(i);
        stats.g2c[i].nq() = stats.g2c[i].tp() + stats.g2c[i].fp();
        stats.g2c[i].fn() = stats.g2c[i].nr() - stats.g2c[i].tp();
    }

    A_ASSERT(stats.oc.nr() >= stats.oc.fn());
 
    for (const auto &i : r.v1())
    {
        if (!stats.findTP(i.name) && isGerm(i))
        {
            VGerm::Match m;
            
            m.var = r.findV1(i.cID, i.l);
            m.rID = i.name;
            A_ASSERT(m.var);
            
            stats.fns.push_back(m);
        }
    }
    
    return stats;
}

static void writeQuins(const FileName &file,
                       const VGerm::SStats &ss,
                       const VGerm::Options &o)
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
    for (const auto &i : r.v1())
    {
        if (isGerm(i))
        {
            // Can we find this sequin?
            const auto isTP = ss.findTP(i.name);
            
            // This shouldn't fail...
            const auto &sv = r.findSeqVar1(i.key());
            
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
                                 % r.af(i.name)
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
                                 % r.af(i.name)
                                 % "-"
                                 % "-"
                                 % gt2str(sv.gt)
                                 % ctx2Str(sv.ctx)
                                 % var2str(i.type())).str());
            }
        }
    }
    
    o.writer->close();
}

static void writeDetected(const FileName &file,
                          const VGerm::SStats &ss,
                          const VGerm::Options &o)
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

    auto f = [&](const std::vector<VGerm::Match> &x, const std::string &label)
    {
        for (const auto &i : x)
        {
            auto sID = (i.var && i.alt && i.ref ? i.var->name : "-");
            const auto ctx = sID != "-" ?  ctx2Str(r.findSeqVar1(i.var->key()).ctx) : "-";

            o.writer->write((boost::format(format) % (i.rID.empty() ? "-" : i.rID)
                                                   % i.qry.cID
                                                   % i.qry.l.start
                                                   % label
                                                   % i.qry.readR
                                                   % i.qry.readV
                                                   % i.qry.depth
                                                   % (sID != "-" ? std::to_string(r.af(sID)) : "-")
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
                         const VGerm::EStats &es,
                         const VGerm::SStats &ss,
                         const VGerm::Options &o)
{
    const auto &r = Standard::instance().r_var;

    extern FileName VCFRef();
    extern FileName BedRef();

    const auto summary = "-------VarGermline Summary Statistics\n\n"
                         "-------VarGermline Output Results\n\n"
                         "       Reference variant annotation:      %1%\n"
                         "       Reference sequin regions:          %2%\n\n"
                         "       Variants identified in sample:     %3%\n"
                         "       Variants identified in sequins:    %4%\n\n"
                         "       Number of sample variants (sequin regions): %5%\n"
                         "       Number of sequin variants (sequin regions): %6%\n\n"
                         "-------Diagnostic performance by variant\n\n"
                         "      *All variants\n"
                         "       Reference:             %7%\n"
                         "       True Positive:         %8%\n"
                         "       False Positive:        %9%\n"
                         "       False Negative:        %10%\n"
                         "       Sensitivity:           %11%\n"
                         "       Precision:             %12%\n"
                         "       F1 Score:              %13%\n"
                         "       FDR Rate:              %14%\n\n"
                         "      *Single Nucleotide Variants (SNVs)\n\n"
                         "       Reference:             %15%\n"
                         "       True Positive:         %16%\n"
                         "       False Positive:        %17%\n"
                         "       False Negative:        %18%\n"
                         "       Sensitivity:           %19%\n"
                         "       Precision:             %20%\n"
                         "       F1 Score:              %21%\n"
                         "       FDR Rate:              %22%\n\n"
                         "      *Small Insertions/Deletions (InDels)\n"
                         "       Reference:             %23%\n"
                         "       True Positive:         %24%\n"
                         "       False Positive:        %25%\n"
                         "       False Negative:        %26%\n"
                         "       Sensitivity:           %27%\n"
                         "       Precision:             %28%\n"
                         "       F1 Score:              %29%\n"
                         "       FDR Rate:              %30%\n\n"
                         "-------Diagnostic performance by genotype\n\n"
                         "      *Homozygous variants\n"
                         "       Reference:             %31%\n"
                         "       True Positive:         %32%\n"
                         "       False Positive:        %33%\n"
                         "       False Negative:        %34%\n"
                         "       Sensitivity:           %35%\n"
                         "       Precision:             %36%\n"
                         "       F1 Score:              %37%\n"
                         "       FDR Rate:              %38%\n\n"
                         "      *Heterozygous variants\n"
                         "       Reference:             %39%\n"
                         "       True Positive:         %40%\n"
                         "       False Positive:        %41%\n"
                         "       False Negative:        %42%\n"
                         "       Sensitivity:           %43%\n"
                         "       Precision:             %44%\n"
                         "       F1 Score:              %45%\n"
                         "       FDR Rate:              %46%\n\n"
                         "-------Diagnostic performance by context\n\n"
                         "       Context                      Sensitivity:\n"
                         "       Common                       %47%\n"
                         "       Low GC                       %48%\n"
                         "       High GC                      %49%\n"
                         "       Long Homopolymer             %50%\n"
                         "       Very Low GC                  %51%\n"
                         "       Very High GC                 %52%\n"
                         "       Short Dinucleotide Repeat    %53%\n"
                         "       Long Dinucleotide Repeat     %54%\n"
                         "       Short Homopolymer            %55%\n"
                         "       Long Quad Nucleotide Repeat  %56%\n"
                         "       Long Trinucleotide Repeat    %57%\n"
                         "       Short Quad Nucleotide Repeat %58%\n"
                         "       Short Trinucleotide Repeat   %59%\n";

        #define D(x) (isnan(x) ? "-" : std::to_string(x))
        
        const auto &snp = ss.v2c.at(Variation::SNP);
        const auto &del = ss.v2c.at(Variation::Deletion);
        const auto &ins = ss.v2c.at(Variation::Insertion);
        const auto &hom = ss.g2c.at(Genotype::Homozygous);
        const auto &het = ss.g2c.at(Genotype::Heterzygous);

        const auto c_nSNP = snp.nq();
        const auto c_nDel = del.nq();
        const auto c_nIns = ins.nq();
        
        auto ind = del;
        ind += ins;

        #define CSN(x) D(ss.c2c.at(x).sn())

        #define E1() (endo.empty() ? "-" : std::to_string(es.v2c.at(Variation::SNP)))
        #define E2() (endo.empty() ? "-" : std::to_string(es.v2c.at(Variation::SNP) + es.v2c.at(Variation::Insertion)))
        #define E3() (endo.empty() ? "-" : std::to_string(es.v2c.at(Variation::SNP) + es.v2c.at(Variation::Insertion) + es.v2c.at(Variation::Deletion)))
        #define E4() (endo.empty() ? "-" : std::to_string(es.g2c.at(Genotype::Homozygous)))
        #define E5() (endo.empty() ? "-" : std::to_string(es.g2c.at(Genotype::Heterzygous)))

        #define C(x) (D(ss.c2c.at(x).nq()))

        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(summary) % VCFRef()                             // 1
                                                % BedRef()                             // 2
                                                % (endo.empty() ? "-" : endo)          // 3
                                                % seqs                                 // 4
                                                % E3()                                 // 5
                                                % (c_nSNP + c_nDel + c_nIns)           // 6
                                                % (r.nType1(Variation::SNP) +
                                                   r.nType1(Variation::Insertion) +
                                                   r.nType1(Variation::Deletion))      // 7
                                                % D(ss.oc.tp())                        // 8
                                                % D(ss.oc.fp())                        // 9
                                                % D(ss.oc.fn())                        // 10
                                                % D(ss.oc.sn())                        // 11
                                                % D(ss.oc.pc())                        // 12
                                                % D(ss.oc.F1())                        // 13
                                                % D(1-ss.oc.pc())                      // 14
                                                % r.nType1(Variation::SNP)             // 15
                                                % D(snp.tp())                          // 16
                                                % D(snp.fp())                          // 17
                                                % D(snp.fn())                          // 18
                                                % D(snp.sn())                          // 19
                                                % D(snp.pc())                          // 20
                                                % D(snp.F1())                          // 21
                                                % D(1 - snp.pc())                      // 22
                                                % (r.nType1(Variation::Insertion) +
                                                   r.nType1(Variation::Deletion))      // 23
                                                % D(ind.tp())                          // 24
                                                % D(ind.fp())                          // 25
                                                % D(ind.fn())                          // 26
                                                % D(ind.sn())                          // 27
                                                % D(ind.pc())                          // 28
                                                % D(ind.F1())                          // 29
                                                % D(1 - ind.pc())                      // 30
                                                % D(r.nGeno1(Genotype::Homozygous))    // 31
                                                % D(hom.tp())                          // 32
                                                % D(hom.fp())                          // 33
                                                % D(hom.fn())                          // 34
                                                % D(hom.sn())                          // 35
                                                % D(hom.pc())                          // 36
                                                % D(hom.F1())                          // 37
                                                % D(1 - hom.pc())                      // 38
                                                % D(r.nGeno1(Genotype::Heterzygous))   // 39
                                                % D(het.tp())                          // 40
                                                % D(het.fp())                          // 41
                                                % D(het.fn())                          // 42
                                                % D(het.sn())                          // 43
                                                % D(het.pc())                          // 44
                                                % D(het.F1())                          // 45
                                                % D(1 - het.pc())                      // 46
                                                % CSN(Context::Common)                 // 47
                                                % CSN(Context::LowGC)                  // 48
                                                % CSN(Context::HighGC)                 // 49
                                                % CSN(Context::LongHompo)              // 50
                                                % CSN(Context::VeryLowGC)              // 51
                                                % CSN(Context::VeryHighGC)             // 52
                                                % CSN(Context::ShortDinRep)            // 53
                                                % CSN(Context::LongDinRep)             // 54
                                                % CSN(Context::ShortHompo)             // 55
                                                % CSN(Context::LongQuadRep)            // 56
                                                % CSN(Context::LongTrinRep)            // 57
                                                % CSN(Context::ShortQuadRep)           // 58
                                                % CSN(Context::ShortTrinRep)           // 59
                         ).str());
    o.writer->close();
}

template <typename T, typename O> void writeFN(const FileName &file, const T &x, const O &o, bool useVar)
{
    extern FileName VCFRef();

    auto wFN = VCFWriter(); wFN.open(o.work + "/" + file);

    std::set<long> keys;
    for (const auto &i : x) { keys.insert(i.var->key()); }

    ParserVCF::parse(VCFRef(), [&](const Variant &x)
    {
        if (keys.count(x.key()))
        {
            wFN.write(x.hdr, x.line);
        }
    });
    
    wFN.close();
}

void VGerm::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto es = analyzeE(endo, o);
    const auto ss = analyzeS(seqs, o);

    o.info("TP: " + std::to_string(ss.oc.tp()));
    o.info("FP: " + std::to_string(ss.oc.fp()));
    o.info("FN: " + std::to_string(ss.oc.fn()));

    o.info("Generating statistics");

    /*
     * Generating VarGermline_sequins.tsv
     */
    
    writeQuins("VarGermline_sequins.tsv", ss, o);

    /*
     * Generating VarGermline_summary.stats
     */
    
    writeSummary("VarGermline_summary.stats", endo, seqs, es, ss, o);
    
    /*
     * Generating VarGermline_detected.tsv
     */
    
    writeDetected("VarGermline_detected.tsv", ss, o);
    
    /*
     * Generating VarGermline_ROC.R
     */
    
    o.generate("VarGermline_ROC.R");
    o.writer->open("VarGermline_ROC.R");
    o.writer->write(createROC("VarGermline_detected.tsv", "data$Depth", "'FP'"));
    o.writer->close();

    /*
     * Generating VarGermline_FN.vcf
     */

    writeFN("VarGermline_FN.vcf", ss.fns, o, true);
}
