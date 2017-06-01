#include "tools/tools.hpp"
#include "VarQuin/v_structure.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_vcf2.hpp"

using namespace Anaquin;

VStructure::EStats VStructure::analyzeE(const FileName &file, const Options &o)
{
    VStructure::EStats stats;
    
    stats.v2c[Variation::SNP];
    stats.v2c[Variation::Deletion];
    stats.v2c[Variation::Inversion];
    stats.v2c[Variation::Insertion];
    stats.v2c[Variation::Duplication];
    
    if (!file.empty())
    {
        ParserVCF2::parse(file, [&](const Variant &x)
        {
            if (o.meth == VStructure::Method::Passed && x.filter != Filter::Pass)
            {
                return;
            }
            else if (!x.opts.count("SVLEN"))
            {
                return;
            }
            
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
        if (o.meth == VStructure::Method::Passed && x.filter != Filter::Pass)
        {
            return;
        }
        else if (!x.opts.count("SVLEN"))
        {
            return;
        }
        
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

    for (const auto &i : r.vars())
    {
        if (!stats.findTP(i.name))
        {
            Match m;
            
            m.var = r.findVar(i.cID, i.l);
            m.rID = i.name;
            A_ASSERT(m.var);
            
            stats.fns.push_back(m);
        }
    }
    
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
    
    const auto summary = "-------VarStructure Output Results\n"
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
                         "       Duplication: %10%\n\n"
                         "-------Sample-derived variants by mutation\n\n"
                         "       Insertion:   %11%\n"
                         "       Deletion:    %12%\n"
                         "       Inversion:   %13%\n"
                         "       Duplication: %14%\n\n"
                         "-------Sequin-derived variants by mutation\n\n"
                         "       Insertion:   %15%\n"
                         "       Deletion:    %16%\n"
                         "       Inversion:   %17%\n"
                         "       Duplication: %18%\n"
                         "       Total:       %19%\n\n"
                         "-------Diagnostic performance by mutation\n\n"
                         "       True Positive:  %20% insertions\n"
                         "       True Positive:  %21% deletions\n"
                         "       True Positive:  %22% inversions\n"
                         "       True Positive:  %23% duplications\n\n"
                         "       False Positive: %24% insertions\n"
                         "       False Positive: %25% deletions\n"
                         "       False Positive: %26% inversions\n"
                         "       False Positive: %27% duplications\n\n"
                         "       False Negative: %28% insertions\n"
                         "       False Negative: %29% deletions\n"
                         "       False Negative: %30% inversions\n"
                         "       False Negative: %31% duplications\n\n"
                         "       *Variants\n"
                         "       Sensitivity: %32$.4f\n"
                         "       Precision:   %33$.4f\n"
                         "       F1 Score:    %34$.4f\n"
                         "       FDR Rate:    %35$.4f\n\n"
                         "       *Insertion\n"
                         "       Sensitivity: %36$.4f\n"
                         "       Precision:   %37$.4f\n"
                         "       F1 Score:    %38$.4f\n"
                         "       FDR Rate:    %39$.4f\n\n"
                         "       *Deletions\n"
                         "       Sensitivity: %40$.4f\n"
                         "       Precision:   %41$.4f\n"
                         "       F1 Score:    %42$.4f\n"
                         "       FDR Rate:    %43$.4f\n\n"
                         "       *Inversions\n"
                         "       Sensitivity: %44$.4f\n"
                         "       Precision:   %45$.4f\n"
                         "       F1 Score:    %46$.4f\n"
                         "       FDR Rate:    %47$.4f\n\n"
                         "       *Duplications\n"
                         "       Sensitivity: %48$.4f\n"
                         "       Precision:   %49$.4f\n"
                         "       F1 Score:    %50$.4f\n"
                         "       FDR Rate:    %51$.4f\n";
    
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
                                                % T()                              // 5
                                                % (nIns + nDel + nInv + nDup)      // 6
                                                % r.nType(Variation::Insertion)    // 7
                                                % r.nType(Variation::Deletion)     // 8
                                                % r.nType(Variation::Inversion)    // 9
                                                % r.nType(Variation::Duplication)  // 10
                                                % E(Variation::Insertion)          // 11
                                                % E(Variation::Deletion)           // 12
                                                % E(Variation::Inversion)          // 13
                                                % E(Variation::Duplication)        // 14
                                                % nIns                             // 15
                                                % nDel                             // 16
                                                % nInv                             // 17
                                                % nDup                             // 18
                                                % (nIns + nDel + nInv + nDup)      // 19
                                                % tp_Ins                           // 20
                                                % tp_Del                           // 21
                                                % tp_Inv                           // 22
                                                % tp_Dup                           // 23
                                                % fp_Ins                           // 24
                                                % fp_Del                           // 25
                                                % fp_Inv                           // 26
                                                % fp_Dup                           // 27
                                                % fn_Ins                           // 28
                                                % fn_Del                           // 29
                                                % fn_Inv                           // 30
                                                % fn_Dup                           // 31
                                                % D(ss.oc.sn())                    // 32
                                                % D(ss.oc.pc())                    // 33
                                                % D(ss.oc.F1())                    // 34
                                                % D(1-ss.oc.pc())                  // 35
                                                % D(ins.sn())                      // 36
                                                % D(ins.pc())                      // 37
                                                % D(ins.F1())                      // 38
                                                % D(1-ins.pc())                    // 39
                                                % D(del.sn())                      // 40
                                                % D(del.pc())                      // 41
                                                % D(del.F1())                      // 42
                                                % D(1-del.pc())                    // 43
                                                % D(inv.sn())                      // 44
                                                % D(inv.pc())                      // 45
                                                % D(inv.F1())                      // 46
                                                % D(1-inv.pc())                    // 47
                                                % D(dup.sn())                      // 48
                                                % D(dup.pc())                      // 49
                                                % D(dup.F1())                      // 50
                                                % D(1-dup.pc())                    // 51
                         ).str());
    o.writer->close();
}

typedef std::vector<Variant> Variants;

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

template <typename T, typename O> void writeVCF(const FileName &file, const T &x, const O &o)
{
    const auto head = "##fileformat=VCFv4.1\n"
                      "##reference=https://www.sequin.xyz\n"
                      "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=""Net Genotype across all datasets"">\n"
                      "##INFO=<ID=END,Number=1,Type=Integer,Description=""End position of the variant described in this record"">\n"
                      "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=""Difference in length between REF and ALT alleles"">\n"
                      "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=""Difference in length between REF and ALT alleles"">\n"
                      "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write(head);

    const auto sv2str = std::map<Variation, std::string>
    {
        { Variation::SNP,         "SNP" },
        { Variation::Deletion,    "DEL" },
        { Variation::Insertion,   "INS" },
        { Variation::Inversion,   "INV" },
        { Variation::Duplication, "DUP" },
    };

    for (const auto &i : x)
    {
        const auto var = i.var ? i.var : &i.qry;
        const auto format = "%1%\t%2%\t%3%\tN\t%4%\t.\t.\tSVTTYPE=%5%;SVLEN=%6%;END=%7%";

        o.writer->write((boost::format(format) % var->cID
                                               % var->l.start
                                               % var->name
                                               % ("<" + sv2str.at(var->type()) + ">")
                                               % sv2str.at(var->type())
                                               % var->opts.at("SVLEN")
                                               % var->l.end).str());
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
    
    writeVCF("VarStructure_TP.vcf", ss.tps, o);

    /*
     * Generating VarStructure_FP.vcf
     */
    
    writeVCF("VarStructure_FP.vcf", ss.fps, o);
    
    /*
     * Generating VarStructure_FN.vcf
     */

    writeVCF("VarStructure_FN.vcf", ss.fns, o);}
