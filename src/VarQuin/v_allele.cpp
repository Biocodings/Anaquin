#include <limits.h>
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_allele.hpp"
#include "parsers/parser_salmon.hpp"

// Defined in resources.cpp
extern Anaquin::FileName MixRef();

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
        
        if (isnan(stats.limit.abund) || exp < stats.limit.abund)
        {
            stats.limit.id = seq.first;
            stats.limit.abund = exp;
        }
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

Scripts VAllele::generateSummary(const FileName &src, const Stats &stats, const VAllele::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto ls = stats.linear();
    
    const auto format = "-------VarKAllele Output\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference VarQuin Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------Allele Frequency\n\n"
                        "       Synthetic: %4%\n"
                        "       Detection Sensitivity: %5% (attomol/ul) (%6%)\n\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %7%\n"
                        "       Correlation: %8%\n"
                        "       R2:          %9%\n"
                        "       F-statistic: %10%\n"
                        "       P-value:     %11%\n"
                        "       SSM:         %12%, DF: %13%\n"
                        "       SSE:         %14%, DF: %15%\n"
                        "       SST:         %16%, DF: %17%\n";
    
    return (boost::format(format) % src                 // 1
                                  % (r.countSeqs()/2.0) // 2
                                  % MixRef()            // 3
                                  % stats.size()        // 4
                                  % stats.limit.abund   // 5
                                  % stats.limit.id      // 6
                                  % ls.m                // 7
                                  % ls.r                // 8
                                  % ls.R2               // 9
                                  % ls.F                // 10
                                  % ls.p                // 11
                                  % ls.SSM              // 12
                                  % ls.SSM_D            // 13
                                  % ls.SSE              // 14
                                  % ls.SSE_D            // 15
                                  % ls.SST              // 16
                                  % ls.SST_D            // 17
                     ).str();
}

void VAllele::writeSummary(const FileName &src, const FileName &file, const Stats &stats, const VAllele::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(VAllele::generateSummary(src, stats, o));
    o.writer->close();
}

void VAllele::report(const FileName &file, const VAllele::Options &o)
{
    const auto stats = analyze(file, o);

    /*
     * Generating VarAllele_summary.csv
     */
    
    VAllele::writeSummary(file, "VarAllele_summary.csv", stats, o);
    
    /*
     * Generating VarAllele_sequins.csv
     */

    VAllele::writeCSV("VarAllele_sequins.csv", stats, o);
    
    /*
     * Generating VarAllele_linear.R
     */

    VAllele::writeRLinear("VarAllele_linear.R", "VarAllele_sequins.csv", stats, o);
}
