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

static void writeQuins(const FileName &file,
                       const VAllele::Stats &stats,
                       const VAllele::Options &o)
{
    const auto format = "%1%\t%2%\t%3%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "ExpFreq"
                                           % "ObsFreq").str());

    const auto d = stats.data(false);
    
    for (auto i = 0; i < d.ids.size(); i++)
    {
        o.writer->write((boost::format(format) % d.ids[i]
                                               % d.x[i]
                                               % d.y[i]).str());
    }
    
    o.writer->close();
}

void VAllele::report(const FileName &file, const VAllele::Options &o)
{
    const auto stats = analyze(file, o);
    
    /*
     * Generating VarAllele_sequins.csv
     */
    
    writeQuins("VarAllele_sequins.csv", stats, o);
    
    /*
     * Generating VarAllele_allele.R
     */
    
    o.generate("VarAllele_allele.R");
    o.writer->open("VarAllele_allele.R");
    o.writer->write(RWriter::createLinear("VarAllele_sequins.csv",
                                          "Allele Frequency",
                                          "Expected allele frequency (log2)",
                                          "Measured allele frequency (log2)",
                                          "ExpFreq",
                                          "ObsFreq",
                                          "input",
                                          true));
    o.writer->close();
}