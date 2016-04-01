#include "VarQuin/v_allele.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotAlleleAllele();

// Defined in resources.cpp
extern Scripts PlotAlleleReads();

static void writeCSV(const FileName &file, const VAllele::Stats &stats, const VAllele::Options &o)
{
    o.writer->open(file);

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";

    auto f = [&](const LinearStats &l, const Label &label)
    {
        const auto data = l.data(false);

        for (auto i = 0; i < data.ids.size(); i++)
        {
            o.writer->write((boost::format(format) % data.ids[i]
                                                   % data.x[i]
                                                   % data.y[i]
                                                   % stats.readR.at(data.ids[i])
                                                   % stats.readV.at(data.ids[i])
                                                   % label).str());
        }
    };

    o.writer->write((boost::format(format) % "Sequin"
                                           % "EAlleleF"
                                           % "MAlleleF"
                                           % "RCount"
                                           % "VCount"
                                           % "Type").str());
    f(stats.snp, "SNP");
    f(stats.ind, "Indel");

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
        if (m.query.cID == ChrT)
        {
            stats.n_chrT++;
            
            switch (m.query.type())
            {
                case Mutation::SNP:       { stats.n_snp++; break; }
                case Mutation::Deletion:
                case Mutation::Insertion: { stats.n_ind++; break; }
            }

            if (m.match && m.ref && m.alt)
            {
                stats.hist.at(m.match->id)++;
                stats.hist.at(toVar(m.match->id))++;

                // Expected allele frequency
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
                stats.all.add(id, known, measured);

                switch (m.query.type())
                {
                    case Mutation::SNP:       { stats.snp.add(id, known, measured); break; }
                    case Mutation::Deletion:
                    case Mutation::Insertion: { stats.ind.add(id, known, measured); break; }
                }

                stats.readR[id] = m.query.readR;
                stats.readV[id] = m.query.readV;
            }
        }
        else
        {
            stats.n_endo++;
        }
    });
    
    stats.all.limit = r.absolute(stats.hist);
    
    return stats;
}

void VAllele::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file, o);
    
    o.info("Detected: " + std::to_string(stats.n_snp) + " SNPs");
    o.info("Detected: " + std::to_string(stats.n_ind) + " indels");

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    o.info("Generating VarAllele_summary.stats");
    o.writer->open("VarAllele_summary.stats");
    o.writer->write(StatsWriter::linearSummary(file, o.rChrT, stats.all, stats.hist));
    o.writer->close();

    /*
     * Generating CSV for all variants
     */

    o.info("Generating VarAllele_quins.csv");
    writeCSV("VarAllele_quins.csv", stats, o);
    
    /*
     * Generating for AlleleReads
     */

    o.info("Generating VarAllele_alleleReads.R");
    o.writer->open("VarAllele_alleleReads.R");
    o.writer->write(RWriter::createScript("VarAllele_quins.csv", PlotAlleleReads()));
    o.writer->close();

    /*
     * Generating for AlleleAllele
     */
    
    o.info("Generating VarAllele_alleleAllele.R");
    o.writer->open("VarAllele_alleleAllele.R");
    o.writer->write(RWriter::createScript("VarAllele_quins.csv", PlotAlleleAllele()));
    o.writer->close();
}