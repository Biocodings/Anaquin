#include <limits.h>
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_allele.hpp"
#include "parsers/parser_salmon.hpp"

using namespace Anaquin;

struct SequinAllele
{
    Measured ref = std::numeric_limits<Measured>::quiet_NaN();
    Measured alt = std::numeric_limits<Measured>::quiet_NaN();
};

VAllele::Stats VAllele::analyze(const FileName &file, const VAllele::Options &o)
{
    VAllele::Stats stats;
    
    std::map<SequinID, SequinAllele> seqs;
    
    switch (o.format)
    {
        case (Format::Salmon):
        {
            ParserSalmon::parse(file, [&](const ParserSalmon::Data &x, const ParserProgress &)
            {
                if (x.cID == ChrIS)
                {
                    if (isRefID(x.id))
                    {
                        seqs[baseID(x.id)].ref = x.abund;
                    }
                    else
                    {
                        seqs[baseID(x.id)].alt = x.abund;
                    }
                }
            });

            break;
        }
    }

    const auto &r = Standard::instance().r_var;

    for (const auto &seq : seqs)
    {
        const auto exp = r.findAFreq(seq.first);
        const auto obs = seq.second.alt / (seq.second.ref + seq.second.alt);
        stats.add(seq.first, exp, obs);
    }

    return stats;
}

Scripts VAllele::generateCSV(const Stats &stats, const VAllele::Options &o)
{
    std::stringstream ss;

    const auto format = "%1%\t%2%\t%3%";
    ss << (boost::format(format) % "ID" % "ExpFreq" % "ObsFreq").str();
    
    const auto d = stats.data(false);
    
    for (auto i = 0; i < d.ids.size(); i++)
    {
        ss << (boost::format(format) % d.ids[i] % d.x[i] % d.y[i]).str() << std::endl;
    }
    
    return ss.str();
}

void VAllele::writeCSV(const FileName &file, const Stats &stats, const VAllele::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(VAllele::generateCSV(stats, o));
    o.writer->close();
}

Scripts VAllele::generateRLinear(const FileName &src, const Stats &stats, const VAllele::Options &o)
{
    return RWriter::createRLinear(src,
                                 o.work,
                                 "Allele Frequency",
                                 "Expected allele frequency (log2)",
                                 "Measured allele frequency (log2)",
                                 "log2(data$ExpFreq)",
                                 "log2(data$ObsFreq)",
                                 "input",
                                 true);
}

void VAllele::writeRLinear(const FileName &file, const FileName &src, const Stats &stats, const VAllele::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(VAllele::generateRLinear(src, stats, o));
    o.writer->close();
}

Scripts VAllele::generateSummary(const Stats &stats, const VAllele::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    
    std::vector<SequinHist>   hists;
    std::vector<LinearStats>  lStats;
    std::vector<MappingStats> mStats;
    
    // Detection limit for the replicates
    Limit limit;
    
    for (auto i = 0; i < files.size(); i++)
    {
        auto &ls = o.metrs == RExpress::Metrics::Isoform ? stats[i].isos : stats[i].genes;
        
        mStats.push_back(stats[i]);
        lStats.push_back(ls);
        
        // Not every replicate is defined...
        if (!stats[i].limit.id.empty())
        {
            if (isnan(limit.abund) || stats[i].limit.abund < limit.abund)
            {
                limit = stats[i].limit;
            }
        }
    }
    
    assert(!isnan(limit.abund) && !limit.id.empty());
    
    const auto title = (o.metrs == Metrics::Gene ? "Genes Expressed" : "Isoform Expressed");
    
    const auto ms = StatsWriter::multiInfect(files, mStats, lStats);
    
    // No reference coordinate annotation given here
    const auto rSyn = o.metrs == Metrics::Gene || shouldAggregate(o) ? r.countGeneSeqs() : r.countSeqs();
    
    const auto format = "-------RnaExpression Output\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference Transcript Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------%4%\n\n"
                        "       Synthetic: %5%\n"
                        "       Detection Sensitivity: %6% (attomol/ul) (%7%)\n\n"
                        "       Genome: %8%\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %9%\n"
                        "       Correlation: %10%\n"
                        "       R2:          %11%\n"
                        "       F-statistic: %12%\n"
                        "       P-value:     %13%\n"
                        "       SSM:         %14%, DF: %15%\n"
                        "       SSE:         %16%, DF: %17%\n"
                        "       SST:         %18%, DF: %19%\n";
    
    return (boost::format(format) % STRING(ms.files)      // 1
                                  % rSyn                  // 2
                                  % MixRef()              // 3
                                  % title                 // 4
                                  % STRING(ms.countSyn)   // 5
                                  % limit.abund           // 6
                                  % limit.id              // 7
                                  % STRING(ms.countGen)   // 8
                                  % STRING(ms.wLog.sl)    // 9
                                  % STRING(ms.wLog.r)     // 10
                                  % STRING(ms.wLog.R2)    // 11
                                  % STRING(ms.wLog.F)     // 12
                                  % STRING(ms.wLog.p)     // 13
                                  % STRING(ms.wLog.SSM)   // 14
                                  % STRING(ms.wLog.SSM_D) // 15
                                  % STRING(ms.wLog.SSE)   // 16
                                  % STRING(ms.wLog.SSE_D) // 17
                                  % STRING(ms.wLog.SST)   // 18
                                  % STRING(ms.wLog.SST_D) // 19
                     ).str();	
}

void VAllele::writeSummary(const FileName &file, const Stats &stats, const VAllele::Options &o)
{
    
}

void VAllele::report(const FileName &file, const VAllele::Options &o)
{
    const auto stats = analyze(file, o);

    /*
     * Generating VarAllele_summary.csv
     */
    
    VAllele::writeSummary("VarAllele_summary.csv", stats, o);
    
    /*
     * Generating VarAllele_sequins.csv
     */

    VAllele::writeCSV("VarAllele_sequins.csv", stats, o);
    
    /*
     * Generating VarAllele_allele.R
     */

    VAllele::writeRLinear("VarAllele_allele.R", "VarAllele_sequins.csv", stats, o);
}
