#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_somatic.hpp"
#include "writers/vcf_writer.hpp"
#include "parsers/parser_vcf.hpp"
#include "writers/json_writer.hpp"

using namespace Anaquin;

typedef SequinVariant::Context Context;

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

static void writeAllele(const VSomatic::Options &o)
{
    extern Scripts PlotAllele();
    
    o.generate("VarMutation_ladder.R");
    o.writer->open("VarMutation_ladder.R");
    o.writer->write(RWriter::createRLinear("VarMutation_sequins.tsv",
                                           o.work,
                                           "Tumor Sample",
                                           "Expected Allele Frequency (log2)",
                                           "Measured Allele Frequency (log2)",
                                           "log2(data$ExpFreq)",
                                           "log2(data$ObsFreq)",
                                           "input",
                                           false,
                                           PlotAllele()));
    o.writer->close();
}

static Scripts createROC(const FileName &file, const std::string &score, const std::string &refRat)
{
    extern Scripts PlotVCROC();
    extern Path __output__;
    extern std::string __full_command__;
    
    return (boost::format(PlotVCROC()) % date()
                                       % __full_command__
                                       % __output__
                                       % file
                                       % score
                                       % refRat).str();
}

EStats VSomatic::analyzeE(const FileName &file, const Options &o)
{
    const auto r2 = Standard::instance().r_var.regs2();
    
    EStats stats;
    
    stats.g2c[Genotype::Homozygous];
    stats.g2c[Genotype::Heterzygous];
    stats.v2c[Variation::SNP];
    stats.v2c[Variation::Deletion];
    stats.v2c[Variation::Inversion];
    stats.v2c[Variation::Insertion];
    stats.v2c[Variation::Duplication];
    
    if (!file.empty())
    {
        auto w = VCFWriter();
        w.open(o.work + "/VarMutation_sample.vcf");

        ParserVCF::parse(file, [&](const Variant &x)
        {
            DInter *d = nullptr;
            
            if (o.filter == VCFFilter::Passed && x.filter != Filter::Pass)
            {
                return;
            }
            else if ((d = contains(r2, x.cID, x.l)))
            {
                w.write(x.hdr, x.line);

                auto t = x;
                t.name = d->name();
                stats.g2c[t.gt]++;
                stats.v2c[t.type()]++;
                stats.vs.insert(t);
            }
        });
        
        w.close();
    }
    
    return stats;
}

/*
 * Strelka implementation
 */

static bool isStrelka(const Variant &x)
{
    auto isSNP = [&]()
    {
        return x.fi.count("AU_2_1") &&
               x.fi.count("CU_2_1") &&
               x.fi.count("GU_2_1") &&
               x.fi.count("TU_2_1") &&
               x.fi.count("AU_2_2") &&
               x.fi.count("CU_2_2") &&
               x.fi.count("GU_2_2") &&
               x.fi.count("TU_2_2");
    };
    
    auto isInd = [&]()
    {
        return x.fi.count("TAR_1_1") &&
               x.fi.count("TAR_1_2") &&
               x.fi.count("TAR_2_1") &&
               x.fi.count("TAR_2_2") &&
               x.fi.count("TIR_1_1") &&
               x.fi.count("TIR_1_2") &&
               x.fi.count("TIR_2_1") &&
               x.fi.count("TIR_2_2");
    };

    return isSNP() || isInd();
}

static Counts strelkaNormalT(const Variant &x, const Sequence &a)
{
    switch (x.type())
    {
        case SNP:
        {
            if (a == "A") { return x.fi.at("AU_1_1"); }
            if (a == "C") { return x.fi.at("CU_1_1"); }
            if (a == "G") { return x.fi.at("GU_1_1"); }
            else          { return x.fi.at("TU_1_1"); }
            break;
        }

        case Deletion:
        case Insertion:
        {
            if (a == x.ref)      { return x.fi.at("TAR_1_1"); }
            else if (a == x.alt) { return x.fi.at("TIR_1_1"); }
            else { throw std::runtime_error("Unknown: " + a); }
        }

        default:        { return NAN; }
    }
}

static Counts strelkaTumorT(const Variant &x, const Sequence &a)
{
    switch (x.type())
    {
        case SNP:
        {
            if (a == "A") { return x.fi.at("AU_2_1"); }
            if (a == "C") { return x.fi.at("CU_2_1"); }
            if (a == "G") { return x.fi.at("GU_2_1"); }
            else          { return x.fi.at("TU_2_1"); }
            break;
        }
            
        case Deletion:
        case Insertion:
        {
            if (a == x.ref)      { return x.fi.at("TAR_2_1"); }
            else if (a == x.alt) { return x.fi.at("TIR_2_1"); }
            else { throw std::runtime_error("Unknown: " + a); }
        }

        default: { return NAN; }
    }
}

// Depth for reference allele in normal sample
static Counts strelkaNormalDPR(const Variant &x)
{
    return strelkaNormalT(x, x.ref);
}

// Depth for variant allele in normal sample
static Counts strelkaNormalDPV(const Variant &x)
{
    return strelkaNormalT(x, x.alt);
}

// Depth for reference allele in tumor sample
static Counts strelkaTumorDPR(const Variant &x)
{
    return strelkaTumorT(x, x.ref);
}

// Depth for variant allele in tumor sample
static Counts strelkaTumorDPV(const Variant &x)
{
    return strelkaTumorT(x, x.alt);
}

static Proportion strelkaNormalAF(const Variant &x)
{
    const auto r = strelkaNormalDPR(x);
    const auto v = strelkaNormalDPV(x);
    return ((Proportion) v) / (r + v);
}

static Proportion strelkaTumorAF(const Variant &x)
{
    const auto r = strelkaTumorDPR(x);
    const auto v = strelkaTumorDPV(x);
    return ((Proportion) v) / (r + v);
}

// Depth for reference allele in normal sample
static Proportion normalDPR(const Variant &x)
{
    if (isStrelka(x)) { return strelkaNormalDPR(x); }
    else              { return NAN; }
}

// Depth for variant allele in normal sample
static Proportion normalDPV(const Variant &x)
{
    if (isStrelka(x)) { return strelkaNormalDPV(x); }
    else              { return NAN; }
}

// Depth for reference allele in tumor sample
static Proportion tumorDPR(const Variant &x)
{
    if (isStrelka(x)) { return strelkaTumorDPR(x); }
    else              { return NAN; }
}

// Depth for variant allele in tumor sample
static Proportion tumorDPV(const Variant &x)
{
    if (isStrelka(x)) { return strelkaTumorDPV(x); }
    else              { return NAN; }
}

// Measured allele frequency for normal
static Proportion normalAF(const Variant &x)
{
    if (x.ff.count("AF_1"))    { return x.ff.at("AF_1");    }
    else if (isStrelka(x))     { return strelkaNormalAF(x); }
    else                       { return NAN; }
}

// Measured allele frequency for tumor
static Proportion tumorAF(const Variant &x)
{
    if (x.ff.count("AF_2"))    { return x.ff.at("AF_2");   }
    else if (x.ff.count("AF")) { return x.ff.at("AF");     }
    else if (isStrelka(x))     { return strelkaTumorAF(x); }
    else                       { return NAN; }
}

VSomatic::SStats VSomatic::analyzeS(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto r2 = Standard::instance().r_var.regs2();
    
    VSomatic::SStats stats;
    
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
        Context::Cancer,
    };
    
    const auto afs = r.vcf1()->lad.groups(Mix_1);
    
    for (auto &i : gts)  { stats.g2c[i]; }
    for (auto &i : ctx)  { stats.c2c[i]; }
    for (auto &i : muts) { stats.v2c[i]; }
    for (auto &i : afs)  { stats.f2c[i]; }
    
    o.analyze(file);
    
    auto wTP = VCFWriter(); wTP.open(o.work + "/VarMutation_TP.vcf");
    auto wFP = VCFWriter(); wFP.open(o.work + "/VarMutation_FP.vcf");
    
    // Caller specific fields
    const auto keys = std::set<std::string>
    {
        "SomaticEVS", "QSS", "QSI",
    };;
    
    ParserVCF::parse(file, [&](const Variant &x)
    {
        if (o.filter == VCFFilter::Passed && x.filter != Filter::Pass)
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
            if (isSomatic(m.var->name))
            {
                for (const auto &i : keys)
                {
                    if (x.fi.count(i))  { stats.si[i][m.var->key()] = x.fi.at(i);  }
                    if (x.ff.count(i))  { stats.sf[i][m.var->key()] = x.ff.at(i);  }
                    if (x.ifi.count(i)) { stats.si[i][m.var->key()] = x.ifi.at(i); }
                    if (x.iff.count(i)) { stats.sf[i][m.var->key()] = x.iff.at(i); }
                }
                
                wTP.write(x.hdr, x.line);
                
                const auto key = m.var->key();
                
                stats.tps.push_back(m);
                
                // Expected allele frequency for tumor
                const auto exp = r.af(m.var->name);
                
                // Measured allele frequency for tumor
                auto obs = tumorAF(x);
                
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
            // Ignore anything for germline
            if (!r.findV2(x.cID, x.l))
            {
                for (const auto &i : keys)
                {
                    if (x.fi.count(i))  { stats.si[i][x.key()] = x.fi.at(i);  }
                    if (x.ff.count(i))  { stats.sf[i][x.key()] = x.ff.at(i);  }
                    if (x.ifi.count(i)) { stats.si[i][x.key()] = x.ifi.at(i); }
                    if (x.iff.count(i)) { stats.sf[i][x.key()] = x.iff.at(i); }
                }
                
                wFP.write(x.hdr, x.line);
                
                // FP because the variant is not found in the reference
                stats.fps.push_back(m);
            }
        }
    });
    
    wTP.close();
    wFP.close();
    
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
            const auto &sv = r.findSeqVar1(i.var->key());
            
            // Overall performance
            stats.oc.tp()++;
            
            // Performance by genotype
            stats.g2c[sv.gt].tp()++;
            
            // Performance by context
            stats.c2c[sv.ctx].tp()++;
            
            // Performance by variation
            stats.v2c[i.qry.type()].tp()++;
            
            // Performance by allele frequency
            stats.f2c[r.af(i.var->name)].tp()++;
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
     * Performance by allele frequency
     */
    
    for (auto &i : r.vcf1()->lad.groups(Mix_1))
    {
        stats.f2c[i].nr() = r.vcf1()->lad.count(i, Mix_1);
        stats.f2c[i].nq() = stats.f2c[i].tp() + stats.f2c[i].fp();
        stats.f2c[i].fn() = stats.f2c[i].nr() - stats.f2c[i].tp();
    }
    
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
        if (!stats.findTP(i.name) && isSomatic(i.name))
        {
            VSomatic::Match m;
            
            m.var = r.findV1(i.cID, i.l);
            m.rID = i.name;
            A_ASSERT(m.var);
            
            stats.fns.push_back(m);
        }
    }
    
    return stats;
}

template <typename T> std::string head(const T &x)
{
    std::stringstream ss;
    
    for (auto &i : x.si) { ss << ("\t" + i.first); }
    for (auto &i : x.sf) { ss << ("\t" + i.first); }
    
    return ss.str();
}

template <typename T> std::string extra(const T &x, long key)
{
    std::stringstream ss;
    
    for (auto &i : x.si)
    {
        if (i.second.count(key)) { ss << ("\t" + toString(i.second.at(key))); }
        else                     { ss << "\t-"; }
    }

    for (auto &i : x.sf)
    {
        if (i.second.count(key)) { ss << ("\t" + toString(i.second.at(key))); }
        else                     { ss << "\t-"; }
    }

    return ss.str();
}

static void writeQuins(const FileName &file,
                       const VSomatic::SStats &ss,
                       const VSomatic::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13$.2f\t%14$.2f\t%15%\t%16%%17%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Position"
                                           % "Label"
                                           % "ReadR_Normal"
                                           % "ReadV_Normal"
                                           % "ReadR_Tumor"
                                           % "ReadV_Tumor"
                                           % "Depth_Normal"
                                           % "Depth_Tumor"
                                           % "ExpFreq"
                                           % "ObsFreq_Normal"
                                           % "ObsFreq_Tumor"
                                           % "Qual"
                                           % "Context"
                                           % "Mutation"
                                           % head(ss)).str());
    for (const auto &i : r.v1())
    {
        if (isSomatic(i.name))
        {
            // Can we find this sequin?
            const auto isTP = ss.findTP(i.name);
            
            // This shouldn't fail...
            const auto &sv = r.findSeqVar1(i.key());

            #define FORMAT_I(x) (c.fi.count(x) ? toString(c.fi.at(x)) : "-")
            #define FORMAT_F(x) (c.ff.count(x) ? toString(c.ff.at(x)) : "-")
            
            if (isTP)
            {
                // Called variant (if found)
                const auto &c = isTP->qry;
                
                o.writer->write((boost::format(format) % i.name
                                                       % i.cID
                                                       % i.l.start
                                                       % "TP"
                                                       % normalDPR(c)
                                                       % normalDPV(c)
                                                       % tumorDPR(c)
                                                       % tumorDPV(c)
                                                       % FORMAT_I("DP_1")
                                                       % FORMAT_I("DP_2")
                                                       % r.af(i.name)
                                                       % normalAF(c)
                                                       % tumorAF(c)
                                                       % toString(c.qual)
                                                       % ctx2Str(sv.ctx)
                                                       % var2str(i.type())
                                                       % extra(ss, i.key())).str());
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
                                                       % "-"
                                                       % "-"
                                                       % "-"
                                                       % r.af(i.name)
                                                       % "-"
                                                       % "-"
                                                       % "-"
                                                       % ctx2Str(sv.ctx)
                                                       % var2str(i.type())
                                                       % extra(ss, i.key())).str());
            }
        }
    }
    
    o.writer->close();
}

static void writeDetected(const FileName &file,
                          const VSomatic::SStats &ss,
                          const VSomatic::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%\t%14%\t%15%\t%16%%17%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Position"
                                           % "Label"
                                           % "ReadR_Normal"
                                           % "ReadV_Normal"
                                           % "ReadR_Somatic"
                                           % "ReadV_Somatic"
                                           % "Depth_Normal"
                                           % "Depth_Somatic"
                                           % "ExpFreq"
                                           % "ObsFreq_Normal"
                                           % "ObsFreq_Tumor"
                                           % "Qual"
                                           % "Context"
                                           % "Mutation"
                                           % head(ss)).str());
    
    auto f = [&](const std::vector<VSomatic::Match> &x, const std::string &label)
    {
        for (const auto &i : x)
        {
            auto sID = (i.var && i.alt && i.ref ? i.var->name : "-");
            const auto ctx = sID != "-" ?  ctx2Str(r.findSeqVar1(i.var->key()).ctx) : "-";
            
            #define _FI_(x) (i.qry.fi.count(x) ? toString(i.qry.fi.at(x)) : "-")
            #define _FF_(x) (i.qry.ff.count(x) ? toString(i.qry.ff.at(x)) : "-")
            
            const auto key = i.var && i.alt && i.ref ? i.var->key() : i.qry.key();
            
            o.writer->write((boost::format(format) % (i.rID.empty() ? "-" : i.rID)
                                                   % i.qry.cID
                                                   % i.qry.l.start
                                                   % label
                                                   % normalDPR(i.qry)
                                                   % normalDPV(i.qry)
                                                   % tumorDPR(i.qry)
                                                   % tumorDPV(i.qry)
                                                   % _FI_("DP_1")
                                                   % _FI_("DP_2")
                                                   % (sID != "-" ? std::to_string(r.af(sID)) : "-")
                                                   % normalAF(i.qry)
                                                   % tumorAF(i.qry)
                                                   % toString(i.qry.qual)
                                                   % ctx
                                                   % var2str(i.qry.type())
                                                   % extra(ss, key)).str());
        }
    };
    
    f(ss.tps, "TP");
    f(ss.fps, "FP");
    
    o.writer->close();
}

static std::map<std::string, std::string> jsonD(const FileName &endo,
                                                const FileName &seqs,
                                                const EStats &es,
                                                const VSomatic::SStats &ss,
                                                const VSomatic::Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    extern bool VCFUser();
    extern FileName VCFRef();
    extern FileName Bed1Ref();

    #define D(x) (isnan(x) ? "-" : std::to_string(x))
    
    const auto &snp = ss.v2c.at(Variation::SNP);
    const auto &del = ss.v2c.at(Variation::Deletion);
    const auto &ins = ss.v2c.at(Variation::Insertion);
    
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
    
    std::map<std::string, std::string> x;
    
    x["vRef"]   = (VCFUser() ? VCFRef() : "-");
    x["bRef"]   = (!Bed1Ref().empty() ? Bed1Ref() : "-");
    x["inputE"] = (endo.empty() ? "-" : endo);
    x["inputS"] = seqs;
    x["nEndo"]  = E3(); // Number of sample variants
    x["nSeqs"]  = D(c_nSNP + c_nDel + c_nIns);
    x["allN"]   = D(r.nType1(Variation::SNP) +
                    r.nType1(Variation::Insertion) +
                    r.nType1(Variation::Deletion));
    x["allTP"]  = D(ss.oc.tp());
    x["allFP"]  = D(ss.oc.fp());
    x["allFN"]  = D(ss.oc.fn());
    x["allSN"]  = D(ss.oc.sn());
    x["allPC"]  = D(ss.oc.pc());
    x["allF1"]  = D(ss.oc.F1());
    x["allFDR"] = D(1-ss.oc.pc());
    x["snpN"]   = D(r.nType1(Variation::SNP));
    x["snpTP"]  = D(snp.tp());
    x["snpFP"]  = D(snp.fp());
    x["snpFN"]  = D(snp.fn());
    x["snpSN"]  = D(snp.sn());
    x["snpPC"]  = D(snp.pc());
    x["snpF1"]  = D(snp.F1());
    x["snpFDR"] = D(1-snp.pc());
    x["indN"]   = D(r.nType1(Variation::Insertion) + r.nType1(Variation::Deletion));
    x["indTP"]  = D(ind.tp());
    x["indFP"]  = D(ind.fp());
    x["indFN"]  = D(ind.fn());
    x["indSN"]  = D(ind.sn());
    x["indPC"]  = D(ind.pc());
    x["indF1"]  = D(ind.F1());
    x["indFDR"] = D(1-ind.pc());
    x["cancer"] = CSN(Context::Cancer);

    return x;
}

static void writeSummary(const FileName &file,
                         const FileName &endo,
                         const FileName &seqs,
                         const EStats &es,
                         const VSomatic::SStats &ss,
                         const VSomatic::Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    extern FileName VCFRef();
    extern FileName Bed1Ref();
    
    std::stringstream str;
    
    const auto grps = r.vcf1()->lad.groups(Mix_1);
    
    for (auto i = grps.rbegin(); i != grps.rend(); i++)
    {
        const auto s = std::to_string(*i);
        str << boost::format("       %1%%2%%3%") % s % std::string(23 - s.size(), ' ') % ss.f2c.at(*i).sn() << std::endl;
    }
    
    const auto summary = "-------VarMutation Summary Statistics\n\n"
                         "-------VarMutation Output Results\n\n"
                         "       Reference variant annotation: %1%\n"
                         "       Reference sequin regions:     %2%\n\n"
                         "       Sample variant file:          %3%\n"
                         "       Sequins variant file:         %4%\n\n"
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
                         "-------Diagnostic performance by context\n\n"
                         "       Context                Sensitivity:\n"
                         "       Cancer                 %31%\n\n"
                         "-------Diagnostic performance by allele frequency\n\n"
                         "       Allele Frequency       Sensitivity:\n" + str.str();

    auto x = jsonD(endo, seqs, es, ss, o);
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % x["vRef"]   // 1
                                            % x["bRef"]   // 2
                                            % x["inputE"] // 3
                                            % x["inputS"] // 4
                                            % x["nEndo"]  // 5
                                            % x["nSeqs"]  // 6
                                            % x["allN"]   // 7
                                            % x["allTP"]  // 8
                                            % x["allFP"]  // 9
                                            % x["allFN"]  // 10
                                            % x["allSN"]  // 11
                                            % x["allPC"]  // 12
                                            % x["allF1"]  // 13
                                            % x["allFDR"] // 14
                                            % x["snpN"]   // 15
                                            % x["snpTP"]  // 16
                                            % x["snpFP"]  // 17
                                            % x["snpFN"]  // 18
                                            % x["snpSN"]  // 19
                                            % x["snpPC"]  // 20
                                            % x["snpF1"]  // 21
                                            % x["snpFDR"] // 22
                                            % x["indN"]   // 23
                                            % x["indTP"]  // 24
                                            % x["indFP"]  // 25
                                            % x["indFN"]  // 26
                                            % x["indSN"]  // 27
                                            % x["indPC"]  // 28
                                            % x["indF1"]  // 29
                                            % x["indFDR"] // 30
                                            % x["cancer"] // 31
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

VSomatic::Stats VSomatic::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    Stats stats;
    
    stats.es = analyzeE(endo, o);
    stats.ss = analyzeS(seqs, o);
    
    o.info("TP: " + std::to_string(stats.ss.oc.tp()));
    o.info("FP: " + std::to_string(stats.ss.oc.fp()));
    o.info("FN: " + std::to_string(stats.ss.oc.fn()));
    
    o.info("Generating statistics");

    /*
     * Generating VarMutation_summary.stats
     */
    
    writeSummary("VarMutation_summary.stats", endo, seqs, stats.es, stats.ss, o);
    
    /*
     * Generating VarMutation_sequins.tsv
     */
    
    writeQuins("VarMutation_sequins.tsv", stats.ss, o);
    
    /*
     * Generating VarMutation_detected.tsv
     */
    
    writeDetected("VarMutation_detected.tsv", stats.ss, o);
    
    /*
     * Generating VarMutation_ROC.R
     */
    
    o.generate("VarMutation_ROC.R");
    o.writer->open("VarMutation_ROC.R");
    o.writer->write(createROC("VarMutation_detected.tsv", "data$ObsFreq_Tumor", "'-'"));
    o.writer->close();
    
    /*
     * Generating VarMutation_ladder.R
     */
    
    writeAllele(o);

    /*
     * Generating VarMutation_FN.vcf
     */
    
    writeFN("VarMutation_FN.vcf", stats.ss.fns, o, true);
    
    /*
     * Generating VarMutation_stats.json
     */
    
    o.generate("VarMutation_stats.json");
    JSONWriter w(o.work);
    w.open("VarMutation_stats.json");
    w.write(jsonD(endo, seqs, stats.es, stats.ss, o));
    w.close();

    return stats;
}
