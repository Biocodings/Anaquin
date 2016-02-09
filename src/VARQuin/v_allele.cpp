#include "VARQuin/v_allele.hpp"

using namespace Anaquin;

static void writeCSV(const FileName &file, const VAllele::Stats &stats, const VAllele::Options &o)
{
    o.writer->open(file);

    const auto format = "%1%\t%2%\t%3%\t%4%";

    auto f = [&](const LinearStats &l, const Label &label)
    {
        const auto data = l.data(false);

        for (auto i = 0; i < data.ids.size(); i++)
        {
            o.writer->write((boost::format(format) % data.ids[i]
                                                   % data.x[i]
                                                   % data.y[i]
                                                   % label).str());
        }
    };

    o.writer->write((boost::format(format) % "Sequin"
                                           % "EAlleleF"
                                           % "MAlleleF"
                                           % "Type").str());
    f(stats.chrT.snp, "SNP");
    f(stats.chrT.ind, "Indel");

    o.writer->close();
}

VAllele::Stats VAllele::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VAllele::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();

    parseVariant(file, o.caller, [&](const VariantMatch &m)
    {
        if (m.query.chrID == ChrT)
        {
            stats.n_chrT++;
            
            if (m.match && m.ref && m.alt)
            {
                stats.hist.at(m.match->id)++;

                // Expected allele frequence
                const auto known = r.alleleFreq(m.match->id);
                
                // Measured coverage is the number of base calls aligned and used in variant calling
                const auto measured = m.query.alleleFreq();

                /*
                 * Plotting the relative allele frequency that is established by differences
                 * in the concentration of reference and variant DNA standards.
                 */
                
                // Eg: D_1_12_R_373892_G/A
                const auto id = (boost::format("%1%_%2%_%3%_%4%:") % m.match->id
                                                                   % m.match->ref
                                                                   % m.match->l.start
                                                                   % m.match->alt).str();
                stats.chrT.tot.add(id, known, measured);

                switch (m.query.type())
                {
                    case Mutation::SNP:       { stats.chrT.snp.add(id, known, measured); break; }
                    case Mutation::Deletion:
                    case Mutation::Insertion: { stats.chrT.ind.add(id, known, measured); break; }
                }
            }
        }
        else
        {
            stats.n_endo++;
            stats.endo.push_back(m.query);
        }
    });
    
    return stats;
}

void VAllele::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    o.writer->open("VarAllele_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT(),
                                                o.rEndo(),
                                                file,
                                                stats,
                                                stats.chrT.tot,
                                                ""));
    o.writer->close();

    /*
     * Generating CSV for all variants
     */

    writeCSV("VarAllele_quins.csv", stats, o);
    
    /*
     * Generating scatter plot for all variants
     */

    o.writer->open("VarAllele_scatter.R");
    o.writer->write(RWriter::scatter(stats.chrT.tot, "Expected vs measured allele frequency", "Expected allele frequency", "Measured allele frequency", "Expected log2 allele frequency", "Measured log2 allele frequency"));
    o.writer->close();

    /*
     * Generating scatter plot for SNPs
     */

    o.writer->open("VarAllele_SNP.R");
    o.writer->write(RWriter::scatter(stats.chrT.snp, "Expected vs measured allele frequency (SNP)", "Expected allele frequency", "Measured allele frequency", "Expected log2 allele frequency", "Measured log2 allele frequency"));
    o.writer->close();
    
    /*
     * Generating scatter plot for indels
     */

    o.writer->open("VarAllele_indel.R");
    o.writer->write(RWriter::scatter(stats.chrT.ind, "Expected vs measured allele frequency (Indel)", "Expected allele frequency", "Measured allele frequency", "Expected log2 allele frequency", "Measured log2 allele frequency"));
    o.writer->close();
}