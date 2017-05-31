#include "tools/tools.hpp"
#include "VarQuin/v_structure.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_vcf2.hpp"

using namespace Anaquin;

extern Scripts PlotVLODR();
extern Scripts PlotVGROC();
extern Scripts PlotVCROC();
extern Scripts PlotAllele();

extern Path __output__;
extern std::string __full_command__;

//static void writeQuins(const FileName &file,
//                       const VStructure::SStats &ss,
//                       const VStructure::Options &o)
//{
//}
//
static void writeDetected(const FileName &file,
                          const VStructure::SStats &ss,
                          const VStructure::Options &o)
{
  
}

static void writeSummary(const FileName &file,
                         const FileName &endo,
                         const FileName &seqs,
                         const VStructure::EStats &es,
                         const VStructure::SStats &ss,
                         const VStructure::Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    extern FileName VCFRef();
    extern FileName BedRef();
    
    const auto summary = "-------VarStructure Output Results\n\n"
                         "-------VarStructure Output\n\n"
                         "       Reference variant annotation:      %1%\n"
                         "       Reference coordinate annotation:   %2%\n\n"
                         "       User identified variants (sample): %3%\n"
                         "       User identified variants (sequin): %4%\n\n"
                         "       Number of sample variants (sequin regions): %5%\n"
                         "       Number of sequin variants (sequin regions): %6%\n\n"
                         "-------Reference variants by mutation\n\n"
                         "       SNPs:   %7%\n"
                         "       Indels: %8%\n"
                         "       Total:  %9%\n\n"
                         "-------Called variants by mutation\n\n"
                         "       %7% SNPs\n"
                         "       %8% indels\n"
                         "       %9% variants\n\n"
                         "-------Diagnostic performance by mutation\n\n"
                         "       True Positive:  %10% Insertion\n"
                         "       True Positive:  %11% Inversion\n"
                         "       True Positive:  %11% Duplication\n"
                         "       True Positive:  %11% Deletion\n\n"
                         "       False Positive:  %10% Insertion\n"
                         "       False Positive:  %11% Inversion\n"
                         "       False Positive:  %11% Duplication\n"
                         "       False Positive:  %11% Deletion\n\n"
                         "       False Negative:  %10% Insertion\n"
                         "       False Negative:  %11% Inversion\n"
                         "       False Negative:  %11% Duplication\n"
                         "       False Negative:  %11% Deletion\n\n"
                         "       *Variants\n"
                         "       Sensitivity: %19$.4f\n"
                         "       Precision:   %20$.4f\n"
                         "       F1 Score:    %21$.4f\n"
                         "       FDR Rate:    %22$.4f\n\n"
                         "       *Insertion\n"
                         "       Sensitivity: %23$.4f\n"
                         "       Precision:   %24$.4f\n"
                         "       F1 Score:    %25$.4f\n"
                         "       FDR Rate:    %26$.4f\n\n"
                         "       *Inversion\n"
                         "       Sensitivity: %27$.4f\n"
                         "       Precision:   %28$.4f\n"
                         "       F1 Score:    %29$.4f\n"
                         "       FDR Rate:    %30$.4f\n\n"
                         "       *Duplication\n"
                         "       Sensitivity: %27$.4f\n"
                         "       Precision:   %28$.4f\n"
                         "       F1 Score:    %29$.4f\n"
                         "       FDR Rate:    %30$.4f\n\n"
                         "       *Deletion\n"
                         "       Sensitivity: %27$.4f\n"
                         "       Precision:   %28$.4f\n"
                         "       F1 Score:    %29$.4f\n"
                         "       FDR Rate:    %30$.4f\n\n";
    
        #define D(x) (isnan(x) ? "-" : std::to_string(x))
        
//        const auto &m2c = ss.m2c;
//        const auto &snp = m2c.at(Variation::SNP);
//        const auto &del = m2c.at(Variation::Deletion);
//        const auto &ins = m2c.at(Variation::Insertion);
//        
//        const auto c_nSNP = snp.nq();
//        const auto c_nDel = del.nq();
//        const auto c_nIns = ins.nq();
//        
//        const auto tp_SNP = snp.tp();
//        const auto tp_Del = del.tp();
//        const auto tp_Ins = ins.tp();
//        
//        const auto fp_SNP = snp.fp();
//        const auto fp_Del = del.fp();
//        const auto fp_Ins = ins.fp();
//        
//        const auto fn_SNP = snp.fn();
//        const auto fn_Del = del.fn();
//        const auto fn_Ins = ins.fn();
//        
//        auto ind = del;
//        ind += ins;
    
        #define CSN(x) D(ss.g2c.at(x).sn())
        
        o.generate(file);
        o.writer->open(file);
//        o.writer->write((boost::format(summary) % VCFRef()                      // 1
//                         % BedRef()                      // 2
//                         % seqs                          // 3
//                         % r.countSNP()                  // 4
//                         % r.countInd()                  // 5
//                         % (r.countSNP() + r.countInd()) // 6
//                         % c_nSNP                        // 7
//                         % (c_nDel + c_nIns)             // 8
//                         % (c_nSNP + c_nDel + c_nIns)    // 9
//                         % tp_SNP                        // 10
//                         % (tp_Del + tp_Ins)             // 11
//                         % (tp_SNP + tp_Del + tp_Ins)    // 12
//                         % fp_SNP                        // 13
//                         % (fp_Del + fp_Ins)             // 14
//                         % (fp_SNP + fp_Del + fp_Ins)    // 15
//                         % fn_SNP                        // 16
//                         % (fn_Del + fn_Ins)             // 17
//                         % (fn_SNP + fn_Del + fn_Ins)    // 18
//                         % D(ss.oc.sn())                 // 19
//                         % D(ss.oc.pc())                 // 20
//                         % D(ss.oc.F1())                 // 21
//                         % D(1-ss.oc.pc())               // 22
//                         % D(snp.sn())                   // 23
//                         % D(snp.pc())                   // 24
//                         % D(snp.F1())                   // 25
//                         % D(1-snp.pc())                 // 26
//                         % D(ind.sn())                   // 27
//                         % D(ind.pc())                   // 28
//                         % D(ind.F1())                   // 29
//                         % D(1-ind.pc())                 // 30
//                         % D(r.nContext(Context::Common))       // 31
//                         % D(r.nContext(Context::VeryLowGC))    // 32
//                         % D(r.nContext(Context::LowGC))        // 33
//                         % D(r.nContext(Context::HighGC))       // 34
//                         % D(r.nContext(Context::VeryHighGC))   // 35
//                         % D(r.nContext(Context::ShortDinRep))  // 36
//                         % D(r.nContext(Context::LongDinRep))   // 37
//                         % D(r.nContext(Context::ShortHompo))   // 38
//                         % D(r.nContext(Context::LongHompo))    // 39
//                         % D(r.nContext(Context::ShortQuadRep)) // 40
//                         % D(r.nContext(Context::LongQuadRep))  // 41
//                         % D(r.nContext(Context::ShortTrinRep)) // 42
//                         % D(r.nContext(Context::LongTrinRep))  // 43
//                         % D(r.nGeno(Genotype::Homozygous))     // 44
//                         % D(r.nGeno(Genotype::Heterzygous))    // 45
//                         % CSN(Context::LowGC)                  // 46
//                         % CSN(Context::HighGC)                 // 47
//                         % CSN(Context::Common)                 // 48
//                         % CSN(Context::LongHompo)              // 49
//                         % CSN(Context::VeryLowGC)              // 50
//                         % CSN(Context::VeryHighGC)             // 51
//                         % CSN(Context::ShortDinRep)            // 52
//                         % CSN(Context::LongDinRep)             // 53
//                         % CSN(Context::ShortHompo)             // 54
//                         % CSN(Context::LongQuadRep)            // 55
//                         % CSN(Context::LongTrinRep)            // 56
//                         % CSN(Context::ShortQuadRep)           // 57
//                         % CSN(Context::ShortTrinRep)           // 58
//                         % (endo.empty() ? "-" : endo)          // 59
//                         % (endo.empty() ? "-" : toString(es.found)) // 60
//                         % (c_nSNP + c_nDel + c_nIns)           // 61
//                         ).str());
    o.writer->close();
}

VStructure::EStats VStructure::analyzeE(const FileName &, const Options &o)
{
    VStructure::EStats stats;
    
    
    return stats;
}

VStructure::SStats VStructure::analyzeS(const FileName &file, const Options &o)
{
    VStructure::SStats stats;
    
    const auto &r = Standard::instance().r_var;
    
    auto vars = std::set<Variation>
    {
        Variation::Deletion,
        Variation::Insertion,
        Variation::Inversion,
        Variation::Duplication,
    };
    
    for (auto &i : vars) { stats.v2c[i]; }
    
    o.analyze(file);

    ParserVCF2::parse(file, [&](const Variant &x)
    {
        const auto &r = Standard::instance().r_var;
        
        auto findMatch = [&](const Variant &qry)
        {
            Match m;
            
            m.qry = qry;
            m.var = nullptr;
            
            // Can we match by position?
            if ((m.var = r.findVar(qry.cID, qry.l)))
            {
                m.rID = m.var->name;
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
        
        if (m.var)
        {
            stats.tps.push_back(m);
        }
        else
        {
            stats.fps.push_back(m);
        }
    });
    
    /*
     * Determining the classification performance
     */
    
    auto forTP = [&]()
    {
        for (auto &i : stats.tps)
        {
            // Overall performance
            stats.oc.tp()++;
            
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
    
    for (auto &var : vars)
    {
        stats.v2c[var].nr() = r.nType(var);
        stats.v2c[var].nq() = stats.v2c[var].tp() + stats.v2c[var].fp();
        stats.v2c[var].fn() = stats.v2c[var].nr() - stats.v2c[var].tp();
        stats.oc.nr() += r.nType(var);
    }

    stats.oc.fn() = stats.oc.nr() - stats.oc.tp();
    
    A_ASSERT(stats.oc.nr() >= stats.oc.fn());
    
    return stats;
}

template <typename Stats, typename Options> void writeQuins(const FileName &file, const Stats &stats, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Start"
                                           % "End"
                                           % "Label"
                                           % "Mutation").str());

    for (const auto &i : r.vars())
    {
        // Can we find this sequin?
        const auto isTP = stats.findTP(i.name);
        
        if (isTP)
        {
            o.writer->write((boost::format(format) % i.name
                                                   % i.cID
                                                   % i.l.start
                                                   % i.l.end
                                                   % "TP"
                                                   % var2str(i.type())).str());
        }
        
        // Failed to detect the variant
        else
        {
            o.writer->write((boost::format(format) % i.name
                                                   % i.cID
                                                   % i.l.start
                                                   % i.l.end
                                                   % "FN"
                                                   % var2str(i.type())).str());
        }
    }
    
    o.writer->close();
}

void VStructure::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto es = analyzeE(endo, o);
    const auto ss = analyzeS(seqs, o);

    /*
     * Generating VarStructure_summary.stats
     */
    
    writeSummary("VarStructure_summary.stats", endo, seqs, es, ss, o);

    /*
     * Generating VarStructure_sequins.csv
     */
    
    writeQuins("VarStructure_sequins.csv", ss, o);
    
    /*
     * Generating VarStructure_detected.csv
     */

    writeDetected("VarStructure_detected.csv", ss, o);
    
    /*
     * Generating VarStructure_TP.vcf
     */

    /*
     * Generating VarStructure_FP.vcf
     */
    
    /*
     * Generating VarStructure_FN.vcf
     */    
}
