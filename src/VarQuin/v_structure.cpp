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
            // This shouldn't fail...
            const auto &sv = r.findSeqVar(i.var->key());
            
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

void VStructure::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto es = analyzeE(endo, o);
    const auto ss = analyzeS(seqs, o);

    /*
     * Generating VarStructure_summary.stats
     */
    
    writeSummary("VarStructure_summary.stats", endo, seqs, es, ss, o);

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
