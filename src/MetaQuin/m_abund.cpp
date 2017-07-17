#include "data/standard.hpp"
#include "MetaQuin/m_abund.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

MAbund::Stats MAbund::analyze(const std::vector<FileName> &files, const MAbund::Options &o)
{
    const auto &r = Standard::instance().r_meta;
 
    MAbund::Stats stats;
    
    for (const auto &i : r.seqsL1())
    {
        stats.hist[i];
    }
    
    switch (o.format)
    {
        case Format::BAM:
        {
            ParserBAM::parse(files[0], [&](const Alignment &x, const ParserBAM::Info &info)
            {
                if (info.p.i && !(info.p.i % 1000000))
                {
                    o.wait(std::to_string(info.p.i));
                }

                if (x.mapped)
                {
                    if (r.seqsL1().count(x.cID))
                    {
                        stats.nSeqs++;
                        stats.hist.at(x.cID)++;
                    }
                    else
                    {
                        stats.nEndo++;
                    }
                }
                else
                {
                    stats.nNA++;
                }
            });

            for (auto &i : stats.hist)
            {
                if (i.second)
                {
                    stats.add(i.first, r.input1(i.first, o.mix), i.second);
                }
            }

            break;
        }

        case Format::RayMeta:
        {
            const auto x = MBlat::analyze(files[1]);
            
            std::map<ContigID, Base> c2l;
            std::map<ContigID, SequinID> c2s;
            std::map<SequinID, std::vector<ContigID>> s2c;
            
            c2l = x.c2l;
            
            for (auto &i : x.aligns)
            {
                c2s[i.first] = i.second->id();
            }
            
            for (auto &i : x.metas)
            {
                for (auto &j : i.second->contigs)
                {
                    s2c[i.first].push_back(j.id);
                }
            }

            // Mapping from contigs to k-mer coverage
            std::map<ContigID, Coverage> c2m;
            
            // Mapping from contigs to k-mer length
            std::map<ContigID, Base> c2kl;
            
            ParserTSV::parse(Reader(files[0]), [&](const ParserTSV::TSV &x)
            {
                c2m[x.id]  = x.kmer;
                c2kl[x.id] = x.klen;
            });
            
            A_ASSERT(!c2m.empty());
            
            /*
             * Quantifying k-mer abundance
             */
            
            for (const auto &i : s2c)
            {
                const auto expected = r.input1(i.first, o.mix);
                auto measured = 0.0;
                
                auto x = 0.0;
                auto y = 0.0;
                
                for (const auto &j : i.second)
                {
                    // K-mer length
                    //const auto kl = c2kl.at(j);
                    
                    // Entire size of the contig (including bases that are not mapped)
                    const auto l = c2l[j];
                    
                    /*
                     * How should we normalize the k-mer observations? Should we normalize by the k-mer length?
                     * Should we normalize by the size of the contig?
                     */
                    
                    x += (double)c2m.at(j); // / kl;// / l;
                    y += l; //j->l.length();
                }
                
                //measured = x / y;
                measured = x;
                
                stats.add(i.first, expected, measured);
            }
            
            break;
        }
    }

    return stats;
}

static void writeQuins(const FileName &file, const MAbund::Stats &stats, const MAbund::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    o.generate(file);
    o.writer->open(file);
    
    std::string format;
    
    switch (o.format)
    {
        case MAbund::Format::BAM:
        {
            format ="%1%\t%2%\t%3%\t%4%\t%5%";
            o.writer->write((boost::format(format) % "ID" % "Length" % "Input" % "Observed" % "FPKM").str());
            break;
        }

        case MAbund::Format::RayMeta:
        {
            format = "%1%\t%2%\t%3%\t%4%";
            o.writer->write((boost::format(format) % "ID" % "Length" % "Input" % "Observed").str());
            break;
        }
    }

    const auto total = sum(stats.hist);
    
    for (const auto &i : stats)
    {
        const auto loc = r.locus(i.first);
        const auto exp = i.second.x;
        const auto obs = i.second.y;
        
        // Normalized measurement
        auto measured = obs;
        
        switch (o.format)
        {
            case MAbund::Format::BAM:
            {
                // Normalized FPKM
                measured = ((double)measured * pow(10, 9)) / (total * loc.length());

                break;
            }

            default: { break; }
        }
        
        switch (o.format)
        {
            case MAbund::Format::BAM:
            {
                o.writer->write((boost::format(format) % i.first
                                                       % loc.length()
                                                       % exp
                                                       % obs
                                                       % measured).str());
                break;
            }

            case MAbund::Format::RayMeta:
            {
                o.writer->write((boost::format(format) % i.first
                                                       % loc.length()
                                                       % exp
                                                       % obs).str());
                break;
            }
        }
    }
    
    o.writer->close();
}

static Scripts generateRLinear(const FileName &src, const MAbund::Stats &stats, const MAbund::Options &o)
{
    switch (o.format)
    {
        case MAbund::Format::BAM:
        {
            return RWriter::createRLinear(src,
                                          o.work,
                                          "FPKM",
                                          "Input Concentration (log2)",
                                          "Measured FPKM (log2)",
                                          "log2(data$Input)",
                                          "log2(data$FPKM)",
                                          "input",
                                          true);
        }

        case MAbund::Format::RayMeta:
        {
            return RWriter::createRLinear(src,
                                          o.work,
                                          "K-mer coverage",
                                          "Input Concentration (log2)",
                                          "Measured K-mer coverage (log2)",
                                          "log2(data$Input)",
                                          "log2(data$Abund)",
                                          "input",
                                          true);
        }
    }
}

static Scripts generateSummary(const FileName &src, const MAbund::Stats &stats, const MAbund::Options &o)
{
    // Defined in resources.cpp
    extern FileName LadRef();

    const auto &r = Standard::instance().r_meta;
    const auto ls = stats.linear();
    
    const auto format = "-------MetaAbund Output\n\n"
                        "       Summary for input: %1%\n\n"
                        "-------Reference MetaQuin Annotations\n\n"
                        "       Synthetic: %2%\n"
                        "       Mixture file: %3%\n\n"
                        "-------Sequin Counts\n\n"
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
    
    const auto limit = stats.limitQuant();
    
    return (boost::format(format) % src               // 1
                                  % r.seqsL1().size() // 2
                                  % LadRef()          // 3
                                  % stats.size()      // 4
                                  % limit.abund       // 5
                                  % limit.id          // 6
                                  % ls.m              // 7
                                  % ls.r              // 8
                                  % ls.R2             // 9
                                  % ls.F              // 10
                                  % ls.p              // 11
                                  % ls.SSM            // 12
                                  % ls.SSM_D          // 13
                                  % ls.SSE            // 14
                                  % ls.SSE_D          // 15
                                  % ls.SST            // 16
                                  % ls.SST_D          // 17
            ).str();
}

static void writeRLinear(const FileName &src, const MAbund::Stats &stats, const MAbund::Options &o)
{
    o.generate("MetaAbund_linear.R");
    o.writer->open("MetaAbund_linear.R");
    o.writer->write(generateRLinear(src, stats, o));
    o.writer->close();
}

void MAbund::report(const std::vector<FileName> &files, const MAbund::Options &o)
{
    const auto stats = MAbund::analyze(files, o);
    
    /*
     * Generating MetaAbund_summary.stats
     */
    
    o.generate("MetaAbund_summary.stats");
    o.writer->open("MetaAbund_summary.stats");
    o.writer->write(generateSummary(files[0], stats, o));
    o.writer->close();
    
    /*
     * Generating MetaAbund_sequins.csv
     */
    
    writeQuins("MetaAbund_sequins.csv", stats, o);
    
    /*
     * Generating MetaAbund_linear.R
     */
    
    writeRLinear("MetaAbund_sequins.csv", stats, o);
}
