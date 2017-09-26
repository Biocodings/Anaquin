#include "VarQuin/VarQuin.hpp"
#include "parsers/parser_vcf.hpp"
#include "VarQuin/v_structure.hpp"

using namespace Anaquin;

VStructure::EStats VStructure::analyzeE(const FileName &file, const Options &o)
{
    VStructure::EStats stats;
    
    stats.v2c[Variation::SNP];
    stats.v2c[Variation::Deletion];
    stats.v2c[Variation::Inversion];
    stats.v2c[Variation::Insertion];
    stats.v2c[Variation::Duplication];
    
    const auto r1 = Standard::instance().r_var.regs1();
    
    if (!file.empty())
    {
        ParserVCF::parse(file, [&](const Variant &x)
        {
            if (o.meth == VStructure::Method::Passed && x.filter != Filter::Pass)
            {
                return;
            }
            else if (!x.ifi.count("SVLEN"))
            {
                return;
            }
            else if (contains(r1, x.cID, x.l))
            {
                stats.v2c[x.type()]++;
            }
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
    
    ParserVCF::parse(file, [&](const Variant &x)
    {
        if (o.meth == VStructure::Method::Passed && x.filter != Filter::Pass)
        {
            return;
        }
        else if (!x.ifi.count("SVLEN"))
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
            if ((m.var = r.findV1(qry.cID, qry.l)))
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
        stats.v2c[var].nr() = r.nType1(var);
        stats.v2c[var].nq() = stats.v2c[var].tp() + stats.v2c[var].fp();
        stats.v2c[var].fn() = stats.v2c[var].nr() - stats.v2c[var].tp();
        stats.oc.nr() += r.nType1(var);
    }
    
    stats.oc.fn() = stats.oc.nr() - stats.oc.tp();
    
    A_ASSERT(stats.oc.nr() >= stats.oc.fn());
    
    for (const auto &i : r.v1())
    {
        if (!stats.findTP(i.name))
        {
            Match m;
            
            m.var = r.findV1(i.cID, i.l);
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
    extern FileName VCFRef();
    extern FileName BedRef();
    
    const auto summary = "-------VarStructure Summary Statistics\n\n"
                         "-------VarStructure Output Results\n\n"
                         "       Reference variant annotation:   %1%\n"
                         "       Reference sequin regions:       %2%\n\n"
                         "       Variants identified in sample:  %3%\n"
                         "       Variants identified in sequins: %4%\n\n"
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
                         "      *Insertions\n"
                         "       Reference:             %15%\n"
                         "       True Positive:         %16%\n"
                         "       False Positive:        %17%\n"
                         "       False Negative:        %18%\n"
                         "       Sensitivity:           %19%\n"
                         "       Precision:             %20%\n"
                         "       F1 Score:              %21%\n"
                         "       FDR Rate:              %22%\n\n"
                         "      *Deletions\n"
                         "       Reference:             %23%\n"
                         "       True Positive:         %24%\n"
                         "       False Positive:        %25%\n"
                         "       False Negative:        %26%\n"
                         "       Sensitivity:           %27%\n"
                         "       Precision:             %28%\n"
                         "       F1 Score:              %29%\n"
                         "       FDR Rate:              %30%\n\n"
                         "      *Inversions\n"
                         "       Reference:             %31%\n"
                         "       True Positive:         %32%\n"
                         "       False Positive:        %33%\n"
                         "       False Negative:        %34%\n"
                         "       Sensitivity:           %35%\n"
                         "       Precision:             %36%\n"
                         "       F1 Score:              %37%\n"
                         "       FDR Rate:              %38%\n\n"
                         "      *Duplications\n"
                         "       Reference:             %39%\n"
                         "       True Positive:         %40%\n"
                         "       False Positive:        %41%\n"
                         "       False Negative:        %42%\n"
                         "       Sensitivity:           %43%\n"
                         "       Precision:             %44%\n"
                         "       F1 Score:              %45%\n"
                         "       FDR Rate:              %46%\n\n";
    
#define D(x) (isnan(x) ? "-" : std::to_string(x))
    
    const auto &ins = ss.v2c.at(Variation::Insertion);
    const auto &del = ss.v2c.at(Variation::Deletion);
    const auto &inv = ss.v2c.at(Variation::Inversion);
    const auto &dup = ss.v2c.at(Variation::Duplication);
    
    const auto nIns = ins.nq();
    const auto nDel = del.nq();
    const auto nInv = inv.nq();
    const auto nDup = dup.nq();
    
#define E(x) (endo.empty() ? "-" : std::to_string(es.v2c.at(x)))
#define T()  (endo.empty() ? "-" : std::to_string(es.v2c.at(Variation::Insertion) + es.v2c.at(Variation::Deletion) + es.v2c.at(Variation::Inversion) + es.v2c.at(Variation::Duplication)))
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % VCFRef()                       // 1
                                            % BedRef()                       // 2
                                            % (endo.empty() ? "-" : endo)    // 3
                                            % seqs                           // 4
                                            % T()                            // 5
                                            % (nIns + nDel + nInv + nDup)    // 6
                                            % D(ss.oc.nr())                  // 7
                                            % D(ss.oc.tp())                  // 8
                                            % D(ss.oc.fp())                  // 9
                                            % D(ss.oc.fn())                  // 10
                                            % D(ss.oc.sn())                  // 11
                                            % D(ss.oc.pc())                  // 12
                                            % D(ss.oc.F1())                  // 13
                                            % D(1-ss.oc.pc())                // 14
                                            % D(ins.nr())                    // 15
                                            % D(ins.tp())                    // 16
                                            % D(ins.fp())                    // 17
                                            % D(ins.fn())                    // 18
                                            % D(ins.sn())                    // 19
                                            % D(ins.pc())                    // 20
                                            % D(ins.F1())                    // 21
                                            % D(1-ins.pc())                  // 22
                                            % D(del.nr())                    // 23
                                            % D(del.tp())                    // 24
                                            % D(del.fp())                    // 25
                                            % D(del.fn())                    // 26
                                            % D(del.sn())                    // 27
                                            % D(del.pc())                    // 28
                                            % D(del.F1())                    // 29
                                            % D(1-del.pc())                  // 30
                                            % D(inv.nr())                    // 31
                                            % D(inv.tp())                    // 32
                                            % D(inv.fp())                    // 33
                                            % D(inv.fn())                    // 34
                                            % D(inv.sn())                    // 35
                                            % D(inv.pc())                    // 36
                                            % D(inv.F1())                    // 37
                                            % D(1-inv.pc())                  // 38
                                            % D(dup.nr())                    // 39
                                            % D(dup.tp())                    // 40
                                            % D(dup.fp())                    // 41
                                            % D(dup.fn())                    // 42
                                            % D(dup.sn())                    // 43
                                            % D(dup.pc())                    // 44
                                            % D(dup.F1())                    // 45
                                            % D(1-dup.pc())                  // 46
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
    
    for (const auto &i : r.v1())
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

template <typename T, typename O> void writeVCF(const FileName &file, const T &x, const O &o, bool useVar)
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
        const auto var = useVar ? i.var : &i.qry;        
        const auto format = "%1%\t%2%\t%3%\tN\t%4%\t%5%\t%6%\t%7%";
        o.writer->write((boost::format(format) % var->cID
                                               % var->l.start
                                               % var->name
                                               % ("<" + sv2str.at(var->type()) + ">")
                                               % "????" //(i.qry.name.empty() ? "." : i.qry.opts.at("QUAL"))
                                               % "????" //(i.qry.name.empty() ? "." : i.qry.opts.at("FILTER"))
                                               % "????" /*(i.qry.name.empty() ? "." : i.qry.opts.at("INFO"))*/).str());
    }
    
    o.writer->close();
}

void VStructure::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto es = analyzeE(endo, o);
    const auto ss = analyzeS(seqs, o);
    
    o.info("TP: " + std::to_string(ss.oc.tp()));
    o.info("FP: " + std::to_string(ss.oc.fp()));
    o.info("FN: " + std::to_string(ss.oc.fn()));
    
    /*
     * Generating VarStructure_summary.stats
     */
    
    writeSummary("VarStructure_summary.stats", endo, seqs, es, ss, o);
    
    /*
     * Generating VarStructure_sequins.csv
     */
    
    writeQuins("VarStructure_sequins.tsv", ss, o);
    
    /*
     * Generating VarStructure_detected.csv
     */
    
    writeDetected("VarStructure_detected.tsv", ss, o);
    
    /*
     * Generating VarStructure_TP.vcf
     */
    
    writeVCF("VarStructure_TP.vcf", ss.tps, o, true);
    
    /*
     * Generating VarStructure_FP.vcf
     */
    
    writeVCF("VarStructure_FP.vcf", ss.fps, o, false);
    
    /*
     * Generating VarStructure_FN.vcf
     */
    
    writeVCF("VarStructure_FN.vcf", ss.fns, o, true);}
