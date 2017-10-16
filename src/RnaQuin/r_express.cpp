#include <cmath>
#include "RnaQuin/r_express.hpp"
#include "parsers/parser_gtf.hpp"
#include "parsers/parser_express.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern FileName LadRef();

struct MultiStats
{
    SStrings files;
    
    SCounts nISeqs, nGSeqs, nIEndo, nGEndo;
    
    // Log-normal regression
    SLinearStats stats;
};

typedef RExpress::Stats   Stats;
typedef RExpress::Format  Format;
typedef RExpress::Options Options;

template <typename F> static void parse(const Reader &r, bool shouldGene, F f)
{
    typedef ParserExpress::Field Field;
    
    ParserProgress p;
    std::string line;
    std::vector<Token> toks;
    
    while (r.nextLine(line))
    {
        Tokens::split(line, "\t", toks);
        ParserExpress::Data x;

        if (p.i)
        {
            x.iID   = toks[Field::IsoID];
            x.cID   = toks[Field::ChrID];
            x.gID   = toks[Field::GeneID];
            x.abund = ss2ld(toks[Field::Abund]);
            f(x, p);
        }
        
        p.i++;
    }
}

template <typename T> void matching(Stats &stats, const T &x, const Options &o)
{
    const auto &r  = Standard::instance().r_rna;
    
    // Sequin transcripts
    const auto &l1 = r.seqsL1();
    
    // Sequin genes
    const auto &l2 = r.seqsL1();

    auto cID = x.cID;
    
    // Fix for quantification programs such as Kallisto
    if (cID == "")
    {
        cID = l1.count(x.iID) || l2.count(x.gID) ? ChrIS() : "endo";
    }

    if (isChrIS(cID))
    {
        SequinID id;
        Concent  exp = NAN;
        Measured obs = NAN;

        Limit *l = nullptr;
        SequinStats *p = nullptr;
        
        // Can we match by isoforms?
        if (l1.count(x.iID))
        {
            if (!isnan(x.abund) && x.abund)
            {
                p   = &stats.isos;
                l   = &stats.iLimit;
                id  = x.iID;
                exp = r.input1(x.iID, o.mix);
                obs = x.abund;
            }
            else
            {
                o.logWarn((boost::format("Zero or invalid measurements for %1%.") % x.iID).str());
            }
            
            if (!isnan(exp)) { stats.nISeqs++; }
        }
        
        // Can we match by genes?
        else if (l2.count(x.gID))
        {
            if (!isnan(x.abund) && x.abund)
            {
                p   = &stats.genes;
                l   = &stats.gLimit;
                id  = x.gID;
                exp = r.input2(x.gID, o.mix);
                obs = x.abund;
            }
            else
            {
                o.logWarn((boost::format("Zero or invalid measurements for %1%.") % x.gID).str());
            }
            
            if (!isnan(exp)) { stats.nGSeqs++; }
        }

        if (p)
        {
            if (p->count(id))
            {
                // This happens to Cufflink guided assembly...
                o.warn("Duplicate: " + id);
                
                p->sum(id, exp, obs);
                
                // Only isoforms are duplicated...
                stats.nISeqs--;
            }
            else
            {
                p->add(id, exp, obs);
            }

            if (isnan(l->abund) || exp < l->abund)
            {
                l->id = id;
                l->abund = exp;
            }
        }
    }
    else
    {
        if (!x.iID.empty()) { stats.nIEndo++; }
        else                { stats.nGEndo++; }

        // We'll need the information to estimate the numbers below and above the LOQ
        stats.gData[x.iID].abund = x.abund;
    }
}

RExpress::Stats RExpress::analyze(const FileName &file, const Options &o)
{
    RExpress::Stats stats;

    o.analyze(file);
    
    switch (o.format)
    {
        case Format::Kallisto:
        {
            ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &x, const ParserProgress &p)
            {
                if (p.i && !(p.i % 100000))
                {
                    o.wait(std::to_string(p.i));
                }
                
                matching(stats, x, o);
            });

            break;
        }
            
        case Format::Anaquin:
        {
            ParserExpress::parse(Reader(file), [&](const ParserExpress::Data &x, const ParserProgress &p)
            {
                if (p.i && !(p.i % 100000))
                {
                    o.wait(std::to_string(p.i));
                }
                
                matching(stats, x, o);
            });

            break;
        }
            
        case Format::GTF:
        {
            ParserExpress::Data t;
            
            ParserGTF::parse(file, [&](const ParserGTF::Data &x, const std::string &, const ParserProgress &p)
            {
                if (p.i && !(p.i % 100000))
                {
                    o.wait(std::to_string(p.i));
                }
                
                t.gID   = x.gID;
                t.iID   = x.tID;
                t.cID   = x.cID;
                t.abund = x.fpkm;

                if (x.type == RNAFeature::Transcript || x.type == RNAFeature::Gene)
                {
                    matching(stats, t, o);
                }
            });

            break;
        }
    }
    
    const auto &r = Standard::instance().r_rna;
    
    // No genes? Only isoforms?
    if (stats.genes.empty() && !stats.isos.empty())
    {
        std::map<GeneID, FPKM> express;
        
        for (const auto &i : stats.isos)
        {
            express[isoform2Gene(i.first)] += i.second.y;
        }
        
        for (const auto &i : express)
        {
            // Input concentration at the gene level
            const auto input = r.input2(i.first, o.mix);
            
            stats.nGSeqs++;
            stats.genes.add(i.first, input, i.second);
            
            if (isnan(stats.gLimit.abund) || input < stats.gLimit.abund)
            {
                stats.gLimit.id = i.first;
                stats.gLimit.abund = input;
            }
        }
    }

    return stats;
}

static Scripts multipleTSV(const std::vector<RExpress::Stats> &stats, bool shouldIso)
{
    const auto &r = Standard::instance().r_rna;
    
    std::set<SequinID> seqs;
    
    // This is the data structure that will be useful
    std::map<unsigned, std::map<SequinID, Concent>> data;
    
    // Expected concentration
    std::map<SequinID, Concent> expect;
    
    std::stringstream ss;
    ss << "ID\tInput\tLength";
    
    for (auto i = 0; i < stats.size(); i++)
    {
        auto &ls = shouldIso ? stats[i].isos : stats[i].genes;

        ss << ((boost::format("\tObserved%1%") % (i+1)).str());

        for (const auto &j : ls)
        {
            seqs.insert(j.first);
            expect[j.first]  = j.second.x;
            data[i][j.first] = j.second.y;
            
            A_CHECK(expect[j.first], "Zero expect concentration");
        }
    }
    
    ss << "\n";
    
    for (const auto &seq : seqs)
    {
        ss << ((boost::format("%1%\t%2%\t%3%") % seq
                                               % expect.at(seq)
                                               % (shouldIso ? r.input3(seq) : r.input4(seq))
                ).str());
        
        for (auto i = 0; i < stats.size(); i++)
        {
            if (data[i].count(seq))
            {
                ss << "\t" << data[i][seq];
            }
            else
            {
                ss << "\tNaN";
            }
        }
        
        ss << "\n";
    }
    
    return ss.str();
}

static MultiStats multiStats(const std::vector<FileName>     &files,
                             const std::vector<RExpress::MappingStats> &mStats,
                             const std::vector<SequinStats>  &lstats)
{
    A_ASSERT(files.size()  == mStats.size());
    A_ASSERT(mStats.size() == lstats.size());

    MultiStats r;
    
    for (auto i = 0; i < lstats.size(); i++)
    {
        const auto lm = lstats[i].linear();
        
        r.nISeqs.add(mStats[i].nISeqs);
        r.nGSeqs.add(mStats[i].nGSeqs);
        r.nIEndo.add(mStats[i].nIEndo);
        r.nGEndo.add(mStats[i].nGEndo);

        r.files.add(files[i]);
        r.stats.p.add(lm.p);
        r.stats.r.add(lm.r);
        r.stats.sl.add(lm.m);
        r.stats.R2.add(lm.R2);
        r.stats.SSM.add(lm.SSM);
        r.stats.SSE.add(lm.SSE);
        r.stats.SST.add(lm.SST);
        r.stats.SSM_D.add(lm.SSM_D);
        r.stats.SSE_D.add(lm.SSE_D);
        r.stats.SST_D.add(lm.SST_D);
    }
    
    return r;
}

Scripts RExpress::generateSummary(const std::vector<FileName> &tmp,
                                  const std::vector<RExpress::Stats> &stats,
                                  const RExpress::Options &o,
                                  const Units &units)
{
    const auto files = path2file(tmp);
    const auto &r = Standard::instance().r_rna;
    
    typedef RExpress::MappingStats MappingStats;
    
    std::vector<SequinStats>  iSStats, gSStats;
    std::vector<MappingStats> iMStats, gMStats;
    
    // Detection limit for the replicates
    Limit iLimit, gLimit;

    auto combine = [&](bool shouldI, std::vector<SequinStats> &ss, std::vector<MappingStats> &ms, Limit &limit)
    {
        for (auto i = 0; i < files.size(); i++)
        {
            auto &ls = shouldI ? stats[i].isos : stats[i].genes;
            
            ss.push_back(ls);
            ms.push_back(stats[i]);
            
            auto p = shouldI ? &stats[i].iLimit : &stats[i].gLimit;
            
            // Not every replicate is defined...
            if (!p->id.empty())
            {
                if (isnan(limit.abund) || p->abund < limit.abund)
                {
                    limit = *p;
                }
            }
        }
    };
    
    combine(true,  iSStats, iMStats, iLimit);
    combine(false, gSStats, gMStats, gLimit);

    A_ASSERT(!isnan(iLimit.abund) && !iLimit.id.empty());
    A_ASSERT(!isnan(gLimit.abund) && !gLimit.id.empty());
    
    const auto iMS = multiStats(files, iMStats, iSStats);
    const auto gMS = multiStats(files, gMStats, gSStats);

    const auto format = "-------RnaExpression Output\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference Transcript Annotations\n\n"
                        "       Sequin: %2% isoforms\n"
                        "       Sequin: %3% genes\n"
                        "       Mixture file: %4%\n\n"
                        "-------Detected Isoforms\n\n"
                        "       Sequin: %5%\n"
                        "       Detection Sensitivity: %6% (attomol/ul) (%7%)\n"
                        "       Genome: %8%\n\n"
                        "-------Linear regression (isoform expression) (log2 scale)\n\n"
                        "       Slope:       %9%\n"
                        "       Correlation: %10%\n"
                        "       R2:          %11%\n"
                        "       P-value:     %12%\n\n"
                        "-------Detected Genes\n\n"
                        "       Sequin: %13%\n"
                        "       Detection Sensitivity: %14% (attomol/ul) (%15%)\n"
                        "       Genome: %16%\n\n"
                        "-------Linear regression (gene expression) (log2 scale)\n\n"
                        "       Slope:       %17%\n"
                        "       Correlation: %18%\n"
                        "       R2:          %19%\n"
                        "       P-value:     %20%\n";
    
    return (boost::format(format) % STRING(iMS.files)    // 1
                                  % r.seqsL1().size()    // 2
                                  % r.seqsL2().size()    // 3
                                  % LadRef()             // 4
                                  % STRING(iMS.nISeqs)   // 5
                                  % iLimit.abund         // 6
                                  % iLimit.id            // 7
                                  % STRING(iMS.nIEndo)   // 8
                                  % STRING(iMS.stats.sl) // 9
                                  % STRING(iMS.stats.r)  // 10
                                  % STRING(iMS.stats.R2) // 11
                                  % STRING(iMS.stats.p)  // 12
                                  % STRING(gMS.nGSeqs)   // 13
                                  % gLimit.abund         // 14
                                  % gLimit.id            // 15
                                  % STRING(gMS.nGEndo)   // 16
                                  % STRING(gMS.stats.sl) // 17
                                  % STRING(gMS.stats.r)  // 18
                                  % STRING(gMS.stats.R2) // 19
                                  % STRING(gMS.stats.p)  // 20
                     ).str();
}

void RExpress::writeSummary(const FileName &file,
                            const std::vector<FileName >&srcs,
                            const std::vector<RExpress::Stats> &stats,
                            const RExpress::Options &o,
                            const Units &units)
{
    o.info("Generating " + file);
    o.writer->open(file);
    o.writer->write(RExpress::generateSummary(srcs, stats, o, units));
    o.writer->close();
}

Scripts RExpress::generateRLinear(const FileName &file,
                                  const std::string &title,
                                  const std::vector<RExpress::Stats> &stats,
                                  const RExpress::Options &o)
{
    Units measured;
    
    switch (o.format)
    {
        case Format::GTF:
        case Format::Anaquin:  { measured = "FPKM"; break; }
        case Format::Kallisto: { measured = "Transcripts per million"; break; }
    }

    if (stats.size() == 1)
    {
        return RWriter::createRLinear(file,
                                      o.work,
                                      title,
                                      "Input Concentration (log2)",
                                       measured + " (log2)",
                                      "log2(data$Input)",
                                      "log2(data$Observed)",
                                      "input",
                                       true);
    }
    else
    {
        return RWriter::createRLinear(file,
                                      o.work,
                                      title,
                                      "Input Concentration (log2)",
                                       measured + " (log2)",
                                      "log2(data$Input)",
                                      "log2(data[,3:ncol(data)])",
                                      "input",
                                       true);
    }
}

void RExpress::writeRLinear(const FileName &file,
                            const FileName &csv,
                            const std::string &title,
                            const std::vector<RExpress::Stats> &stats,
                            const RExpress::Options &o)
{
    o.info("Generating " + file);
    o.writer->open(file);
    o.writer->write(RExpress::generateRLinear(csv, title, stats, o));
    o.writer->close();
}

Scripts RExpress::generateITSV(const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    const auto &r = Standard::instance().r_rna;
    
    if (stats.size() == 1)
    {
        std::stringstream ss;
        
        const auto format = "%1%\t%2%\t%3%\t%4%\n";
        ss << (boost::format(format) % "Name" % "Length" % "Input" % "Observed").str();
        
        for (const auto &i : stats[0].isos)
        {
            ss << (boost::format(format) % i.first
                                         % r.input3(i.first)
                                         % i.second.x
                                         % i.second.y).str();
        }

        return ss.str();
    }
    else
    {
        return multipleTSV(stats, true);
    }
}

Scripts RExpress::generateGTSV(const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    const auto &r = Standard::instance().r_rna;
    
    if (stats.size() == 1)
    {
        std::stringstream ss;
        
        const auto format = "%1%\t%2%\t%3%\t%4%\n";
        ss << (boost::format(format) % "Name" % "Length" % "Input" % "Observed").str();
        
        for (const auto &i : stats[0].genes)
        {
            ss << (boost::format(format) % i.first
                                         % r.input4(i.first)
                                         % i.second.x
                                         % i.second.y).str();
        }
        
        return ss.str();
    }
    else
    {
        return multipleTSV(stats, false);
    }
}

void RExpress::writeITSV(const FileName &output, const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    o.info("Generating " + output);
    o.writer->open(output);
    o.writer->write(generateITSV(stats, o));
    o.writer->close();
}

void RExpress::writeGTSV(const FileName &output, const std::vector<RExpress::Stats> &stats, const RExpress::Options &o)
{
    o.info("Generating " + output);
    o.writer->open(output);
    o.writer->write(generateGTSV(stats, o));
    o.writer->close();
}

void RExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = analyze(files, o);

    /*
     * Generating RnaExpression_summary.stats
     */
    
    RExpress::writeSummary("RnaExpression_summary.stats", files, stats, o, "isoforms");
    
    /*
     * Generating RnaExpression_isoforms.tsv
     */
    
    RExpress::writeITSV("RnaExpression_isoforms.tsv", stats, o);

    /*
     * Generating RnaExpression_genes.tsv
     */

    RExpress::writeGTSV("RnaExpression_genes.tsv", stats, o);
    
    /*
     * Generating RnaExpression_isoforms.R
     */
    
    RExpress::writeRLinear("RnaExpression_isoforms.R", "RnaExpression_isoforms.tsv", "Isoform Expression", stats, o);

    /*
     * Generating RnaExpression_genes.R
     */

    RExpress::writeRLinear("RnaExpression_genes.R", "RnaExpression_genes.tsv", "Gene Expression", stats, o);
}
