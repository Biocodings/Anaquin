#include "VarQuin/v_allele.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVAllele();

// Defined in resources.cpp
extern Scripts PlotVAlleleReads();

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

    o.writer->write((boost::format(format) % "sequin"
                                           % "expected"
                                           % "measured"
                                           % "rcount"
                                           % "vcount"
                                           % "type").str());
    
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

    switch (o.soft)
    {
        case Software::Kallisto:
        {
            /*
             * The implementation differs to a variant caller. Typically, we'd estimate by
             * the number of reads supporting the reference and alternative allele. Obviously,
             * we don't have the alleles here. We should model by pooling the reference and
             * variant sequins.
             */
            
            std::set<SequinID> ids;
            std::map<SequinID, Coverage> matchr, matchv;

            ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &d, const ParserProgress &)
            {
                const auto m = r.match(d.id);
                
                if (m)
                {
                    const auto bID = baseID(d.id);
                    ids.insert(bID);
                    
                    if (isRefID(d.id))
                    {
                        matchr[bID] = d.abund;
                    }
                    else
                    {
                        matchv[bID] = d.abund;
                    }
                    
                    stats.n_chrT++;
                    stats.hist.at(m->id)++;
                }
            });
            
            for (const auto &id : ids)
            {
                if (matchr.count(id) && matchv.count(id))
                {
                    const auto ref = matchr[id];
                    const auto var = matchv[id];
                    
                    // Expected abundance
                    const auto known = r.matchAlleleFreq(id);
                    
                    // Measured abundance
                    const auto measured = var / (ref + var);
                    
                    stats.all.add(id, known, measured);
                }
            }
            
            assert(stats.snp.empty());
            assert(stats.ind.empty());
            assert(stats.readR.empty());
            assert(stats.readV.empty());
            
            break;
        }

        default:
        {
            parseVariant(file, o.soft, [&](const VariantMatch &m)
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
                        
                        // Expected allele frequency
                        const auto known = r.matchAlleleFreq(baseID(m.match->id));
                        
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
                    stats.n_geno++;
                }
            });

            break;
        }
    }

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
    o.writer->write(StatsWriter::linearSummary(file, o.rChrT, stats.all, stats, stats.hist, "variants"));
    o.writer->close();

    /*
     * Generating detailed statistics
     */

    o.info("Generating VarAllele_quins.csv");
    
    switch (o.soft)
    {
        case VAllele::Software::Kallisto:
        {
            o.writer->open("VarAllele_quins.csv");
            o.writer->write(StatsWriter::writeCSV(stats.all));
            o.writer->close();
            break;
        }

        default:
        {
            writeCSV("VarAllele_quins.csv", stats, o);

            /*
             * Generating for allele vs reads
             */
            
            o.info("Generating VarAllele_reads.R");
            o.writer->open("VarAllele_reads.R");
            o.writer->write(RWriter::createScript("VarAllele_quins.csv", PlotVAlleleReads()));
            o.writer->close();
            
            break;
        }
    }

    /*
     * Generating for allele vs allele
     */
    
    o.info("Generating VarAllele_allele.R");
    o.writer->open("VarAllele_allele.R");
    o.writer->write(RWriter::createScript("VarAllele_quins.csv", PlotVAllele()));
    o.writer->close();
    
    /*
     * Generating a report
     */
    
    o.report->open("VAllele_report.pdf");
    o.report->addTitle("VAllele_report");
    o.report->addFile("VarAllele_summary.stats");
    o.report->addFile("VarAllele_quins.csv");
    o.report->addFile("VarAllele_allele.R");
    o.report->addFile("VarAllele_reads.R");
}