#include "VarQuin/v_freq.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

static void writeQuins(const FileName &file, const VFreq::Stats &stats, const VFreq::Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";
    const auto x = stats.vars.data(false);

    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Position"
                                           % "Expected"
                                           % "Observed"
                                           % "ReadR"
                                           % "ReadV"
                                           % "Depth"
                                           % "Type").str());
    
    for (const auto &i : stats.hist)
    {
        const auto &cID = i.first;
        
        if (Standard::isSynthetic(cID))
        {
            for (const auto &j : i.second)
            {
                const auto m = r.findVar(cID, j.first);
                assert(m);

                const auto &key = toString(j.first);
                
                const auto type  = type2str(m->type());
                const auto uname = (m->id + "_" + toString(m->l.start) + "_" + type);
                
                if (x.id2x.count(key))
                {
                    o.writer->write((boost::format(format) % uname
                                                           % m->l.start
                                                           % x.id2x.at(key)
                                                           % x.id2y.at(key)
                                                           % stats.readR.at(j.first)
                                                           % stats.readV.at(j.first)
                                                           % stats.depth.at(j.first)
                                                           % type).str());
                }
                else
                {
                    o.writer->write((boost::format(format) % uname
                                                           % m->l.start
                                                           % "NA"
                                                           % "NA"
                                                           % "NA"
                                                           % "NA"
                                                           % "NA"
                                                           % type).str());
                }
            }
        }
    }

    o.writer->close();
}

VFreq::Stats VFreq::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VFreq::Stats stats;
    
    stats.hist = r.vHist();
    
    for (const auto &i : stats.hist)
    {
        stats.data[i.first];
    }
    
    o.info("Reading VCF inputs");
    stats.vData = vcfData(file, o.input);
    
    o.info("Parsing VCF inputs");

    parseVariants(file, o.input, [&](const VariantMatch &m)
    {
        const auto &cID = m.query.cID;

        auto f = [&]()
        {
            // Only matching if the position and alleles agree
            const auto matched = m.match && m.ref && m.alt_;
            
            if (matched)
            {
                auto key = m.match->key(); //var2hash(m.match->id, m.match->type(), m.match->l);
                stats.hist.at(cID).at(key)++;
                //stats.hist.at(m.match->id)++;
                
                stats.data.at(cID).af = m.query.alleleFreq();
                
                if (Standard::isSynthetic(cID))
                {
                    const auto exp = r.findAFreq(baseID(m.match->id));
                    const auto obs = m.query.alleleFreq();
                 
                    // Eg: 2821292107
                    const auto id = toString(key);
                    
                    // Add for all variants
                    stats.vars.add(id, exp, obs);
                    
                    switch (m.query.type())
                    {
                        case Mutation::SNP:       { stats.snp.add(id, exp, obs); break; }
                        case Mutation::Deletion:
                        case Mutation::Insertion: { stats.ind.add(id, exp, obs); break; }
                    }
                    
                    stats.readR[key] = m.query.readR;
                    stats.readV[key] = m.query.readV;
                    stats.depth[key] = m.query.depth;
                    
                    if (isnan(stats.vars.limit.abund) || exp < stats.vars.limit.abund)
                    {
                        stats.vars.limit.id = m.match->id;
                        stats.vars.limit.abund = exp;
                    }
                }
            }
        };
 
        if (Standard::isSynthetic(cID))
        {
            stats.n_syn++;
            f();
        }
        else
        {
            stats.n_gen++;
            
            if (Standard::isGenomic(cID))
            {
                f();
            }
        }
    });

    return stats;
}

static Scripts generateSummary(const FileName &file, const VFreq::Stats &stats, const VFreq::Options &o)
{
    const auto &r = Standard::instance().r_var;

    extern FileName VCFRef();
    extern FileName MixRef();

    const auto lm = stats.vars.linear(true);
    
    // Calcluate the inflection point with logarithm
    auto ms = stats.vars.limitQuant(true);
    
    // Remember the break-point is on the log2-scale, we'll need to convert it back
    ms.b = pow(2, ms.b);
    
    Counts n_below = 0;
    Counts n_above = 0;
    
    for (const auto &i : stats.data)
    {
        if (!Standard::isSynthetic(i.first))
        {
            if (i.second.af >= ms.b)
            {
                n_above++;
            }
            else
            {
                n_below++;
            }
        }
    }

    const auto mm = r.findVar(ChrT, stol(ms.id));
    assert(mm);
    
    const auto hasGen = n_below || n_above;

    #define G(x) (hasGen ? toString(x) : "-")

    const auto summary = "-------VarFrequency Output\n\n"
                         "      Reference variant annotations: %1%\n"
                         "      User identified variants: %2%\n"
                         "      Sequin mixture file: %3%\n\n"
                         "-------Reference variant annotations\n\n"
                         "      Synthetic: %4% SNPs\n"
                         "      Synthetic: %5% indels\n"
                         "      Synthetic: %6% variants\n\n"
                         "      Genome: %7% SNPs\n"
                         "      Genome: %8% indels\n"
                         "      Genome: %9% variants\n\n"
                         "-------User variant identification\n\n"
                         "      Synthetic: %10% SNPs\n"
                         "      Synthetic: %11% indels\n"
                         "      Synthetic: %12% variants\n\n"
                         "      Detection Sensitivity: %13% (attomol/ul) (%14%)\n\n"
                         "      Genome: %15% SNPs\n"
                         "      Genome: %16% indels\n"
                         "      Genome: %17% variants\n\n"
                         "-------Limit of Quantification (LOQ)\n"
                         "      *Estimated by piecewise segmented regression\n\n"
                         "       Break: %18% attomol/ul (%19%)\n\n"
                         "      *Below LOQ\n"
                         "       Intercept:   %20%\n"
                         "       Slope:       %21%\n"
                         "       Correlation: %22%\n"
                         "       R2:          %23%\n"
                         "       Genome:      %24%\n\n"
                         "      *Above LOQ\n"
                         "       Intercept:   %25%\n"
                         "       Slope:       %26%\n"
                         "       Correlation: %27%\n"
                         "       R2:          %28%\n"
                         "       Genome:      %29%\n\n"
                         "-------Overall linear regression (log2 scale)\n\n"
                         "      Correlation: %30%\n"
                         "      Slope:       %31%\n"
                         "      R2:          %32%\n"
                         "      F-statistic: %33%\n"
                         "      P-value:     %34%\n"
                         "      SSM:         %35%, DF: %36%\n"
                         "      SSE:         %37%, DF: %38%\n"
                         "      SST:         %39%, DF: %40%\n";
    
    return ((boost::format(summary) % file                         // 1
                                    % VCFRef()                     // 2
                                    % MixRef()                     // 3
                                    % r.countSNPSyn()              // 4
                                    % r.countIndSyn()              // 5
                                    % r.countVarSyn()              // 6
                                    % r.countSNPGen()              // 7
                                    % r.countIndGen()              // 8
                                    % r.countVarGen()              // 9
                                    % stats.vData.countSNPSyn()    // 10
                                    % stats.vData.countIndSyn()    // 11
                                    % stats.vData.countVarSyn()    // 12
                                    % stats.vars.limit.abund       // 13
                                    % stats.vars.limit.id          // 14
                                    % G(stats.vData.countSNPGen()) // 15
                                    % G(stats.vData.countIndGen()) // 16
                                    % G(stats.vData.countVarGen()) // 17
                                    % ms.b                         // 18
                                    % mm->id                       // 19
                                    % ms.lInt                      // 20
                                    % ms.lSl                       // 21
                                    % ms.lr                        // 22
                                    % ms.lR2                       // 23
                                    % n_above                      // 24
                                    % ms.rInt                      // 25
                                    % ms.rSl                       // 26
                                    % ms.rr                        // 27
                                    % ms.rR2                       // 28
                                    % n_below                      // 29
                                    % lm.r                         // 30
                                    % lm.m                         // 31
                                    % lm.R2                        // 32
                                    % lm.F                         // 33
                                    % lm.p                         // 34
                                    % lm.SSM                       // 35
                                    % lm.SSM_D                     // 36
                                    % lm.SSE                       // 37
                                    % lm.SSE_D                     // 38
                                    % lm.SST                       // 39
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
     * Generating VarFrequency_sequins.csv
     */

    o.info("Generating VarFrequency_sequins.csv");
    writeQuins("VarFrequency_sequins.csv", stats, o);
    
    /*
     * Generating VarFrequency_allele.R
     */
    
    o.info("Generating VarFrequency_allele.R");
    o.writer->open("VarFrequency_allele.R");
    o.writer->write(RWriter::createScatterNeedLog("VarFrequency_sequins.csv",
                                                  "Allele Frequency",
                                                  "Expected allele frequency (log2)",
                                                  "Measured allele frequency (log2)",
                                                  "Expected",
                                                  "Observed", true));
    o.writer->close();
    
    /*
     * Generating VarFrequency_reads.R
     */
    
    //o.info("Generating VarFrequency_reads.R");
    //o.writer->open("VarFrequency_reads.R");
    //o.writer->write(RWriter::createScript("VarFrequency_sequins.csv", PlotScatter()));
    //o.writer->close();
    
    /*
     * Generating VarFrequency_report.pdf
     */
    
    o.report->open("VarFrequency_report.pdf");
    o.report->addTitle("VarFrequency_report");
    o.report->addFile("VarFrequency_summary.stats");
    o.report->addFile("VarFrequency_sequins.csv");
    o.report->addFile("VarFrequency_allele.R");
    o.report->addFile("VarFrequency_reads.R");
}