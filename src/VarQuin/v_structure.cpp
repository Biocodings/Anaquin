#include "tools/tools.hpp"
#include "VarQuin/v_structure.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_vcf2.hpp"

using namespace Anaquin;

VStructure::EStats VStructure::analyzeE(const FileName &file, const Options &)
{
    VStructure::EStats stats;
    
    stats.v2c[Variation::Insertion];
    stats.v2c[Variation::Deletion];
    stats.v2c[Variation::Inversion];
    stats.v2c[Variation::Duplication];
    
    if (!file.empty())
    {
        ParserVCF2::parse(file, [&](const Variant &x)
        {
            stats.v2c[x.type()]++;
        });
    }

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

static void writeDetected(const FileName &file,
                          const VStructure::SStats &ss,
                          const VStructure::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Position"
                                           % "Label"
                                           % "Mutation").str());
    
    auto f = [&](const std::vector<VStructure::Match> &x, const std::string &label)
    {
        for (const auto &i : x)
        {
            o.writer->write((boost::format(format) % (i.rID.empty() ? "-" : i.rID)
                                                   % i.qry.cID
                                                   % i.qry.l.start
                                                   % label
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
                         "       Insertion:   %7%\n"
                         "       Deletion:    %8%\n"
                         "       Inversion:   %9%\n"
                         "       Duplication: %10%\n"
                         "       Total:       %11%\n\n"
                         "-------Sample-derived variants by mutation\n\n"
                         "       Insertion:   %12%\n"
                         "       Deletion:    %13%\n"
                         "       Inversion:   %14%\n"
                         "       Duplication: %15%\n"
                         "       Total:       %16%\n\n"
                         "-------Sequin-derived variants by mutation\n\n"
                         "       Insertion:   %17%\n"
                         "       Deletion:    %18%\n"
                         "       Inversion:   %19%\n"
                         "       Duplication: %20%\n"
                         "       Total:       %21%\n\n"
                         "-------Diagnostic performance by mutation\n\n"
                         "       True Positive:  %22% insertions\n"
                         "       True Positive:  %23% deletions\n"
                         "       True Positive:  %24% inversions\n"
                         "       True Positive:  %25% duplications\n\n"
                         "       False Positive: %26% insertions\n"
                         "       False Positive: %27% deletions\n"
                         "       False Positive: %28% inversions\n"
                         "       False Positive: %29% duplications\n\n"
                         "       False Negative: %30% insertions\n"
                         "       False Negative: %31% deletions\n"
                         "       False Negative: %32% inversions\n"
                         "       False Negative: %33% duplications\n\n"
                         "       *Variants\n"
                         "       Sensitivity: %34$.4f\n"
                         "       Precision:   %35$.4f\n"
                         "       F1 Score:    %36$.4f\n"
                         "       FDR Rate:    %37$.4f\n\n"
                         "       *Insertion\n"
                         "       Sensitivity: %38$.4f\n"
                         "       Precision:   %39$.4f\n"
                         "       F1 Score:    %40$.4f\n"
                         "       FDR Rate:    %41$.4f\n\n"
                         "       *Deletions\n"
                         "       Sensitivity: %42$.4f\n"
                         "       Precision:   %43$.4f\n"
                         "       F1 Score:    %44$.4f\n"
                         "       FDR Rate:    %45$.4f\n\n"
                         "       *Inversions\n"
                         "       Sensitivity: %46$.4f\n"
                         "       Precision:   %47$.4f\n"
                         "       F1 Score:    %48$.4f\n"
                         "       FDR Rate:    %49$.4f\n\n"
                         "       *Duplications\n"
                         "       Sensitivity: %50$.4f\n"
                         "       Precision:   %51$.4f\n"
                         "       F1 Score:    %52$.4f\n"
                         "       FDR Rate:    %53$.4f\n";
    
        #define D(x) (isnan(x) ? "-" : std::to_string(x))
        
        const auto &ins = ss.v2c.at(Variation::Insertion);
        const auto &del = ss.v2c.at(Variation::Deletion);
        const auto &inv = ss.v2c.at(Variation::Inversion);
        const auto &dup = ss.v2c.at(Variation::Duplication);
    
        const auto nIns = ins.nq();
        const auto nDel = del.nq();
        const auto nInv = inv.nq();
        const auto nDup = dup.nq();

        const auto tp_Ins = ins.tp();
        const auto tp_Del = del.tp();
        const auto tp_Inv = inv.tp();
        const auto tp_Dup = dup.tp();
    
        const auto fp_Ins = ins.fp();
        const auto fp_Del = del.fp();
        const auto fp_Inv = inv.fp();
        const auto fp_Dup = dup.fp();
    
        const auto fn_Ins = ins.fn();
        const auto fn_Del = del.fn();
        const auto fn_Inv = inv.fn();
        const auto fn_Dup = dup.fn();
    
        #define E(x) (endo.empty() ? "-" : std::to_string(es.v2c.at(x)))
        #define T()  (endo.empty() ? "-" : std::to_string(es.v2c.at(Variation::Insertion) + es.v2c.at(Variation::Deletion) + es.v2c.at(Variation::Inversion) + es.v2c.at(Variation::Duplication)))

        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(summary) % VCFRef()                         // 1
                                                % BedRef()                         // 2
                                                % (endo.empty() ? "-" : endo)      // 3
                                                % seqs                             // 4
                                                % (endo.empty() ? "-" : endo)      // 5
                                                % (nIns + nDel + nInv + nDup)      // 6
                                                % r.nType(Variation::Insertion)    // 7
                                                % r.nType(Variation::Deletion)     // 8
                                                % r.nType(Variation::Inversion)    // 9
                                                % r.nType(Variation::Duplication)  // 10
                                                % (r.nType(Variation::Insertion) +
                                                   r.nType(Variation::Deletion)  +
                                                   r.nType(Variation::Inversion) +
                                                   r.nType(Variation::Duplication))
                                                % E(Variation::Insertion)           // 12
                                                % E(Variation::Deletion)            // 13
                                                % E(Variation::Inversion)           // 14
                                                % E(Variation::Duplication)         // 15
                                                % T()                               // 16
                                                % nIns                              // 17
                                                % nDel                              // 18
                                                % nInv                              // 19
                                                % nDup                              // 20
                                                % (nIns + nDel + nInv + nDup)       // 21
                                                % tp_Ins                            // 22
                                                % tp_Del                            // 23
                                                % tp_Inv                            // 24
                                                % tp_Dup                            // 25
                                                % fp_Ins                            // 26
                                                % fp_Del                            // 27
                                                % fp_Inv                            // 28
                                                % fp_Dup                            // 29
                                                % fn_Ins                            // 30
                                                % fn_Del                            // 31
                                                % fn_Inv                            // 32
                                                % fn_Dup                            // 33
                                                % D(ss.oc.sn())                     // 34
                                                % D(ss.oc.pc())                     // 35
                                                % D(ss.oc.F1())                     // 36
                                                % D(1-ss.oc.pc())                   // 37
                                                % D(ins.sn())                       // 38
                                                % D(ins.pc())                       // 39
                                                % D(ins.F1())                       // 40
                                                % D(1-ins.pc())                     // 41
                                                % D(del.sn())                       // 42
                                                % D(del.pc())                       // 43
                                                % D(del.F1())                       // 44
                                                % D(1-del.pc())                     // 45
                                                % D(inv.sn())                       // 46
                                                % D(inv.pc())                       // 47
                                                % D(inv.F1())                       // 48
                                                % D(1-inv.pc())                     // 49
                                                % D(dup.sn())                       // 50
                                                % D(dup.pc())                       // 51
                                                % D(dup.F1())                       // 52
                                                % D(1-dup.pc())                     // 53
                         ).str());
    o.writer->close();
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
