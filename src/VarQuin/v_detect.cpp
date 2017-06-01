#include "tools/tools.hpp"
#include "VarQuin/v_detect.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_vcf2.hpp"

using namespace Anaquin;

typedef SequinVariant::Context Context;

extern Scripts PlotVGROC();
extern Path __output__;
extern std::string __full_command__;

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

    stats.g2c[Genotype::Homozygous];
    stats.g2c[Genotype::Heterzygous];
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
            
            stats.g2c[x.gt]++;
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

    for (auto &i : gts)  { stats.g2c[i]; }
    for (auto &i : ctx)  { stats.c2c[i]; }
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
        stats.v2c[mut].nr() = r.nType(mut);
        stats.v2c[mut].nq() = stats.v2c[mut].tp() + stats.v2c[mut].fp();
        stats.v2c[mut].fn() = stats.v2c[mut].nr() - stats.v2c[mut].tp();
        stats.oc.nr() += r.nType(mut);
    }

    stats.oc.fn() = stats.oc.nr() - stats.oc.tp();
    
    /*
     * Performance by context
     */
    
    for (auto &i : ctx)
    {
        stats.c2c[i].nr() = r.nContext(i);
        stats.c2c[i].nq() = stats.c2c[i].tp() + stats.c2c[i].fp();
        stats.c2c[i].fn() = stats.c2c[i].nr() - stats.c2c[i].tp();
    }
    
    /*
     * Performance by genotype
     */
    
    for (auto &i : gts)
    {
        stats.g2c[i].nr() = r.nGeno(i);
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

    const auto summary = "-------VarDetect Output Results\n\n"
                         "-------VarDetect Output\n\n"
                         "       Reference variant annotation:      %1%\n"
                         "       Reference coordinate annotation:   %2%\n\n"
                         "       User identified variants (sample): %3%\n"
                         "       User identified variants (sequin): %4%\n\n"
                         "       Number of reference variants:      %5%\n"
                         "       Number of sample-derived variants: %6%\n"
                         "       Number of sequin-derived variants: %7%\n\n"
                         "-------Reference variants by mutation\n\n"
                         "       SNPs:   %8%\n"
                         "       Indels: %9%\n\n"
                         "-------Reference variants by genotype\n\n"
                         "       Homozygosity:   %10%\n"
                         "       Heterozygosity: %11%\n\n"
                         "-------Reference variants by context\n\n"
                         "       Common:                       %12%\n"
                         "       Very Low GC:                  %13%\n"
                         "       Low GC:                       %14%\n"
                         "       High GC:                      %15%\n"
                         "       Very High GC:                 %16%\n"
                         "       Short Dinucleotide Repeat:    %17%\n"
                         "       Long Dinucleotide Repeat:     %18%\n"
                         "       Short Homopolymer:            %19%\n"
                         "       Long Homopolymer:             %20%\n"
                         "       Short Quad Nucleotide Repeat: %21%\n"
                         "       Long Quad Nucleotide Repeat:  %22%\n"
                         "       Short Trinucleotide Repeat:   %23%\n"
                         "       Long Trinucleotide Repeat:    %24%\n\n"
                         "-------Sample-derived variants by mutation\n\n"
                         "       %25% SNPs\n"
                         "       %26% indels\n\n"
                         "-------Sequin-derived variants by mutation\n\n"
                         "       %27% SNPs\n"
                         "       %28% indels\n\n"
                         "-------Sample-derived variants by genotype\n\n"
                         "       %29% Homozygosity\n"
                         "       %30% Heterozygosity\n\n"
                         "-------Sequin-derived variants by genotype\n\n"
                         "       %31% Homozygosity\n"
                         "       %32% Heterozygosity\n\n"
                         "-------Sequin-derived variants by context\n\n"
                         "       %33% Low GC\n"
                         "       %34% High GC\n"
                         "       %35% Common\n"
                         "       %36% Long Homopolymer\n"
                         "       %37% Very Low GC\n"
                         "       %38% Very High GC\n"
                         "       %39% Short Dinucleotide Repeat\n"
                         "       %40% Long Dinucleotide Repeat\n"
                         "       %41% Short Homopolymer\n"
                         "       %42% Long Quad Nucleotide Repeat\n"
                         "       %43% Long Trinucleotide Repeat\n"
                         "       %44% Short Quad Nucleotide Repeat\n"
                         "       %45% Short Trinucleotide Repeat\n\n"
                         "-------Overall diagnostic performance\n\n"
                         "       Sensitivity: %46$.4f\n"
                         "       Precision:   %47$.4f\n"
                         "       F1 Score:    %48$.4f\n"
                         "       FDR Rate:    %49$.4f\n\n"
                         "-------Diagnostic performance by mutation\n\n"
                         "       True Positive:  %50% SNPs\n"
                         "       True Positive:  %51% indels\n"
                         "       True Positive:  %52% variants\n\n"
                         "       False Positive: %53% SNPs\n"
                         "       False Positive: %54% indels\n"
                         "       False Positive: %55% variants\n\n"
                         "       False Negative: %56% SNPs\n"
                         "       False Negative: %57% indels\n"
                         "       False Negative: %58% variants\n\n"
                         "       *SNPs\n"
                         "       Sensitivity: %59$.4f\n"
                         "       Precision:   %60$.4f\n"
                         "       F1 Score:    %61$.4f\n"
                         "       FDR Rate:    %62$.4f\n\n"
                         "       *Indels\n"
                         "       Sensitivity: %63$.4f\n"
                         "       Precision:   %64$.4f\n"
                         "       F1 Score:    %65$.4f\n"
                         "       FDR Rate:    %66$.4f\n\n"
                         "-------Diagnostic performance by genotype\n\n"
                         "       True Positive:  %67% Homozygosity\n"
                         "       True Positive:  %68% Heterozygosity\n\n"
                         "       False Positive: %69% Homozygosity\n"
                         "       False Positive: %70% Heterozygosity\n\n"
                         "       False Negative: %71% Homozygosity\n"
                         "       False Negative: %72% Heterozygosity\n\n"
                         "       *Homozygosity\n"
                         "       Sensitivity: %73$.4f\n"
                         "       Precision:   %74$.4f\n"
                         "       F1 Score:    %75$.4f\n"
                         "       FDR Rate:    %76$.4f\n\n"
                         "       *Heterozygosity\n"
                         "       Sensitivity: %77$.4f\n"
                         "       Precision:   %78$.4f\n"
                         "       F1 Score:    %79$.4f\n"
                         "       FDR Rate:    %80$.4f\n\n"
                         "-------Diagnostic performance by context\n\n"
                         "       *Low GC\n"
                         "       Sensitivity: %81$.4f\n\n"
                         "       *High GC\n"
                         "       Sensitivity: %82$.4f\n\n"
                         "       *Common\n"
                         "       Sensitivity: %83$.4f\n\n"
                         "       *Long Homopolymer\n"
                         "       Sensitivity: %84$.4f\n\n"
                         "       *Very Low GC\n"
                         "       Sensitivity: %85$.4f\n\n"
                         "       *Very High GC\n"
                         "       Sensitivity: %86$.4f\n\n"
                         "       *Short Dinucleotide Repeat\n"
                         "       Sensitivity: %87$.4f\n\n"
                         "       *Long Dinucleotide Repeat\n"
                         "       Sensitivity: %88$.4f\n\n"
                         "       *Short Homopolymer\n"
                         "       Sensitivity: %89$.4f\n\n"
                         "       *Long Quad Nucleotide Repeat\n"
                         "       Sensitivity: %90$.4f\n\n"
                         "       *Long Trinucleotide Repeat\n"
                         "       Sensitivity: %91$.4f\n\n"
                         "       *Short Quad Nucleotide Repeat\n"
                         "       Sensitivity: %92$.4f\n\n"
                         "       *Short Trinucleotide Repeat\n"
                         "       Sensitivity: %93$.4f";

        #define D(x) (isnan(x) ? "-" : std::to_string(x))
        
        const auto &snp = ss.v2c.at(Variation::SNP);
        const auto &del = ss.v2c.at(Variation::Deletion);
        const auto &ins = ss.v2c.at(Variation::Insertion);
        const auto &hom = ss.g2c.at(Genotype::Homozygous);
        const auto &het = ss.g2c.at(Genotype::Heterzygous);

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

        #define CSN(x) D(ss.c2c.at(x).sn())

        #define E1() (endo.empty() ? "-" : std::to_string(es.v2c.at(Variation::SNP)))
        #define E2() (endo.empty() ? "-" : std::to_string(es.v2c.at(Variation::SNP) + es.v2c.at(Variation::Insertion)))
        #define E3() (endo.empty() ? "-" : std::to_string(es.v2c.at(Variation::SNP) + es.v2c.at(Variation::Insertion) + es.v2c.at(Variation::Deletion)))
        #define E4() (endo.empty() ? "-" : std::to_string(es.g2c.at(Genotype::Homozygous)))
        #define E5() (endo.empty() ? "-" : std::to_string(es.g2c.at(Genotype::Heterzygous)))

        #define C(x) (D(ss.c2c.at(Context::LowGC).nq()))

        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(summary) % VCFRef()                             // 1
                                                % BedRef()                             // 2
                                                % seqs                                 // 3
                                                % (endo.empty() ? "-" : endo)          // 4
                                                % (r.nType(Variation::SNP) +
                                                   r.nType(Variation::Insertion) +
                                                   r.nType(Variation::Deletion))       // 5
                                                % E3()                                 // 6
                                                % (c_nSNP + c_nDel + c_nIns)           // 7
                                                % r.nType(Variation::SNP)              // 8
                                                % (r.nType(Variation::Insertion) +
                                                   r.nType(Variation::Deletion))       // 9
                                                % D(r.nGeno(Genotype::Homozygous))     // 10
                                                % D(r.nGeno(Genotype::Heterzygous))    // 11
                                                % D(r.nContext(Context::Common))       // 12
                                                % D(r.nContext(Context::VeryLowGC))    // 13
                                                % D(r.nContext(Context::LowGC))        // 14
                                                % D(r.nContext(Context::HighGC))       // 15
                                                % D(r.nContext(Context::VeryHighGC))   // 16
                                                % D(r.nContext(Context::ShortDinRep))  // 17
                                                % D(r.nContext(Context::LongDinRep))   // 18
                                                % D(r.nContext(Context::ShortHompo))   // 19
                                                % D(r.nContext(Context::LongHompo))    // 20
                                                % D(r.nContext(Context::ShortQuadRep)) // 21
                                                % D(r.nContext(Context::LongQuadRep))  // 22
                                                % D(r.nContext(Context::ShortTrinRep)) // 23
                                                % D(r.nContext(Context::LongTrinRep))  // 24
                                                % E1()                                 // 25
                                                % E2()                                 // 26
                                                % c_nSNP                               // 27
                                                % (c_nDel + c_nIns)                    // 28
                                                % E4()                                 // 29
                                                % E5()                                 // 30
                                                % D(ss.g2c.at(Genotype::Homozygous).nq())                               // 31
                                                % D(ss.g2c.at(Genotype::Heterzygous).nq())                               // 32
                                                % (C(Context::LowGC))                  // 33
                                                % (C(Context::HighGC))                 // 34
                                                % (C(Context::Common))                 // 35
                                                % (C(Context::LongHompo))              // 36
                                                % (C(Context::VeryLowGC))              // 37
                                                % (C(Context::VeryHighGC))             // 38
                                                % (C(Context::ShortDinRep))            // 39
                                                % (C(Context::LongDinRep))             // 40
                                                % (C(Context::ShortHompo))             // 41
                                                % (C(Context::LongHomp))               // 42
                                                % (C(Context::LongQuadRep))            // 43
                                                % (C(Context::ShortTrinRep))           // 44
                                                % (C(Context::LongTrinRep))            // 45
                                                % D(ss.oc.sn())                        // 46
                                                % D(ss.oc.pc())                        // 47
                                                % D(ss.oc.F1())                        // 48
                                                % D(1-ss.oc.pc())                      // 49
                                                % tp_SNP                               // 50
                                                % (tp_Del + tp_Ins)                    // 51
                                                % (tp_SNP + tp_Del + tp_Ins)           // 52
                                                % fp_SNP                               // 53
                                                % (fp_Del + fp_Ins)                    // 54
                                                % (fp_SNP + fp_Del + fp_Ins)           // 55
                                                % fn_SNP                               // 56
                                                % (fn_Del + fn_Ins)                    // 57
                                                % (fn_SNP + fn_Del + fn_Ins)           // 58
                                                % D(snp.sn())                          // 59
                                                % D(snp.pc())                          // 60
                                                % D(snp.F1())                          // 61
                                                % D(1-snp.pc())                        // 62
                                                % D(ind.sn())                          // 63
                                                % D(ind.pc())                          // 64
                                                % D(ind.F1())                          // 65
                                                % D(1-ind.pc())                        // 66
                                                % D(hom.tp())                          // 67
                                                % D(het.tp())                          // 68
                                                % D(hom.fp())                          // 69
                                                % D(het.fp())                          // 70
                                                % D(hom.fn())                          // 71
                                                % D(het.fn())                          // 72
                                                % D(hom.sn())                          // 73
                                                % D(hom.pc())                          // 74
                                                % D(hom.F1())                          // 75
                                                % D(1-hom.pc())                        // 76
                                                % D(het.sn())                          // 77
                                                % D(het.pc())                          // 78
                                                % D(het.F1())                          // 79
                                                % D(1-het.pc())                        // 80
                                                % CSN(Context::LowGC)                  // 81
                                                % CSN(Context::HighGC)                 // 82
                                                % CSN(Context::Common)                 // 83
                                                % CSN(Context::LongHompo)              // 84
                                                % CSN(Context::VeryLowGC)              // 85
                                                % CSN(Context::VeryHighGC)             // 86
                                                % CSN(Context::ShortDinRep)            // 87
                                                % CSN(Context::LongDinRep)             // 88
                                                % CSN(Context::ShortHompo)             // 89
                                                % CSN(Context::LongQuadRep)            // 90
                                                % CSN(Context::LongTrinRep)            // 91
                                                % CSN(Context::ShortQuadRep)           // 92
                                                % CSN(Context::ShortTrinRep)           // 93
                         ).str());
    
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
        const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t.\t.\t.";
        
        o.writer->write((boost::format(format) % var->cID
                                           % var->l.start
                                           % var->name
                                           % var->ref
                                           % var->alt).str());
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
     * Generating VarDa _TP.vcf
     */
    
    writeVCF("VarDetect_TP.vcf", ss.tps, o);
    
    /*
     * Generating VarDetect_FP.vcf
     */
    
    writeVCF("VarDetect_FP.vcf", ss.fps, o);
    
    /*
     * Generating VarDetect_FN.vcf
     */
    
    writeVCF("VarDetect_FN.vcf", ss.fns, o);
}
