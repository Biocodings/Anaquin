#include "VarQuin/v_freq.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

static void writeCSV(const FileName &file, const VFreq::Stats &stats, const VFreq::Options &o)
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

    o.writer->write((boost::format(format) % "Seq"
                                           % "Expected"
                                           % "Observed"
                                           % "ReadsR"
                                           % "ReadsV"
                                           % "Type").str());
    
    f(stats.snp, "SNP");
    f(stats.ind, "Indel");

    o.writer->close();
}

VFreq::Stats VFreq::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VFreq::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();

    switch (o.input)
    {
//        case Input::Kallisto:
//        {
//            /*
//             * The implementation differs to a variant caller. Typically, we'd estimate by
//             * the number of reads supporting the reference and alternative allele. Obviously,
//             * we don't have the alleles here. We should model by pooling the reference and
//             * variant sequins.
//             */
//            
//            std::set<SequinID> ids;
//            std::map<SequinID, Coverage> matchr, matchv;
//
//            ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &d, const ParserProgress &)
//            {
//                const auto m = r.match(d.id);
//                
//                if (m)
//                {
//                    const auto bID = baseID(d.id);
//                    ids.insert(bID);
//                    
//                    if (isRefID(d.id))
//                    {
//                        matchr[bID] = d.abund;
//                    }
//                    else
//                    {
//                        matchv[bID] = d.abund;
//                    }
//                    
//                    stats.n_syn++;
//                    stats.hist.at(m->id)++;
//                }
//            });
//            
//            for (const auto &id : ids)
//            {
//                if (matchr.count(id) && matchv.count(id))
//                {
//                    const auto ref = matchr[id];
//                    const auto var = matchv[id];
//                    
//                    // Expected abundance
//                    const auto known = r.findAFreq(id);
//                    
//                    // Measured abundance
//                    const auto measured = var / (ref + var);
//                    
//                    stats.all.add(id, known, measured);
//                }
//            }
//            
//            assert(stats.snp.empty());
//            assert(stats.ind.empty());
//            assert(stats.readR.empty());
//            assert(stats.readV.empty());
//            
//            break;
//        }

        default:
        {
            parseVariants(file, o.input, [&](const VariantMatch &m)
            {
                if (m.query.cID == ChrT)
                {
                    stats.n_syn++;
                    
                    switch (m.query.type())
                    {
                        case Mutation::SNP:       { stats.n_snp++; break; }
                        case Mutation::Deletion:
                        case Mutation::Insertion: { stats.n_ind++; break; }
                    }
                    
                    if (m.match && m.ref && m.alt)
                    {
                        stats.hist.at(m.match->id)++;
                        
                        const auto expected = r.findAFreq(baseID(m.match->id));
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
                        stats.vars.add(id, expected, measured);

                        switch (m.query.type())
                        {
                            case Mutation::SNP:       { stats.snp.add(id, expected, measured); break; }
                            case Mutation::Deletion:
                            case Mutation::Insertion: { stats.ind.add(id, expected, measured); break; }
                        }

                        stats.readR[id] = m.query.readR;
                        stats.readV[id] = m.query.readV;
                    }
                }
                else
                {
                    stats.n_gen++;
                }
            });

            break;
        }
    }

    stats.vars.limit = r.absolute(stats.hist);

    return stats;
}

static Scripts generateSummary(const FileName &file, const VFreq::Stats &stats, const VFreq::Options &o)
{
    const auto &r = Standard::instance().r_var;

    extern FileName VCFRef();
    extern FileName MixRef();

    const auto lm = stats.vars.linear(true);

    const auto summary = "-------VarFrequency Output\n\n"
                         "Reference variant annotations: %1%\n"
                         "User identified variants: %2%\n"
                         "Sequin mixture file: %3%\n\n"
                         "-------Reference variant annotations\n\n"
                         "Synthetic: %4% SNPs\n"
                         "Synthetic: %5% indels\n"
                         "Synthetic: %6% variants\n\n"
                         "-------User Variant Identification\n\n"
                         "Synthetic: %7%\n"
                         "Detection Sensitivity: %8% (attomol/ul) (%9%)\n\n"
                         "-------Overall linear regression (log2 scale)\n\n"
                         "Correlation: %10%\n"
                         "Slope:       %11%\n"
                         "R2:          %12%\n"
                         "F-statistic: %13%\n"
                         "P-value:     %14%\n"
                         "SSM:         %15%, DF: %16%\n"
                         "SSE:         %17%, DF: %18%\n"
                         "SST:         %19%, DF: %20%\n";
    
    return ((boost::format(summary) % file
                                    % VCFRef()
                                    % MixRef()
                                    % r.countSNPSync()
                                    % r.countIndSync()
                                    % r.countSync()
                                    % "????" // 7
                                    % stats.vars.limit.abund // 8
                                    % stats.vars.limit.id    // 9
                                    % lm.r // 10
                                    % lm.m
                                    % lm.R2
                                    % lm.F
                                    % lm.p
                                    % lm.SSM
                                    % lm.SSM_D
                                    % lm.SSE
                                    % lm.SSE_D
                                    % lm.SST
                                    % lm.SST_D).str());
}

void VFreq::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file, o);
    
    o.info("Generating statistics");

    /*
     * Generating VarFrequency_summary.stats
     */

    o.info("Generating VarFrequency_summary.stats");
    o.writer->open("VarFrequency_summary.stats");
    o.writer->write(generateSummary(file, stats, o));
    o.writer->close();

    /*
     * Generating VarFrequency_quins.csv
     */

    o.info("Generating VarFrequency_quins.csv");
    writeCSV("VarFrequency_quins.csv", stats, o);
    
    /*
     * Generating VarFrequency_allele.R
     */
    
    o.info("Generating VarFrequency_allele.R");
    o.writer->open("VarFrequency_allele.R");
    o.writer->write(RWriter::createScatter("VarFrequency_quins.csv",
                                           "Allele Frequency",
                                           "Expected allele frequency (log2)",
                                           "Measured allele frequency (log2)",
                                           "Expected",
                                           "Observed"));
    o.writer->close();
    
    /*
     * Generating VarFrequency_reads.R
     */
    
    //o.info("Generating VarFrequency_reads.R");
    //o.writer->open("VarFrequency_reads.R");
    //o.writer->write(RWriter::createScript("VarFrequency_quins.csv", PlotScatter()));
    //o.writer->close();
    
    /*
     * Generating VarFrequency_report.pdf
     */
    
    o.report->open("VarFrequency_report.pdf");
    o.report->addTitle("VarFrequency_report");
    o.report->addFile("VarFrequency_summary.stats");
    o.report->addFile("VarFrequency_quins.csv");
    o.report->addFile("VarFrequency_allele.R");
    o.report->addFile("VarFrequency_reads.R");
}