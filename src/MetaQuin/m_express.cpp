#include "parsers/parser_sam.hpp"
#include "MetaQuin/m_express.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMExpress();

MExpress::Stats MExpress::analyze(const std::vector<FileName> &files, const MExpress::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MExpress::Stats stats;
    stats.hist = r.hist();
    
    for (auto &file : files)
    {
        o.analyze(file);
        
        switch (o.soft)
        {
            case MExpress::BWA:
            case MExpress::Bowtie:
            {
                ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
                {
                    if (align.cID == ChrT)
                    {
                        throw std::runtime_error("Invalid alignments. Please check the documentation and try again.");
                    }
                    
                    if (align.mapped)
                    {
                        const auto m = r.match(align.cID);
                        
                        if (m)
                        {
                            stats.hist.at(m->id)++;
                        }
                    }
                });
                
                break;
            }
                
            case MExpress::Velvet:
            {
                break;
            }
                
            case MExpress::RayMeta:
            {
                /*
                 *      Format: <Contigs> <TSV> <PSL>
                 */

                const auto t = MBlat::analyze(files[2]);

                /*
                 * Mapping contigs to k-mer coverage
                 */
                
                std::map<ContigID, Coverage> c2m;

                ParserTSV::parse(Reader(files[1]), [&](const ParserTSV::TSV &x, const ParserProgress &)
                {
                    c2m[x.id] = x.kmer;
                });
                
                assert(!c2m.empty());

                /*
                 * Mapping sequins to k-mer coverage
                 */
                
                std::map<ContigID, SequinID> c2s;
                std::map<SequinID, std::vector<const AlignedContig *>> s2c;

                for (auto &seq : stats.hist)
                {
                    if (t.metas.count(seq.first))
                    {
                        for (auto &i : t.metas.at(seq.first)->contigs)
                        {
                            s2c[seq.first].push_back(&i);
                            c2s[i.id] = seq.first;
                        }
                    }
                }
                
                /*
                 * Quantifying expressions
                 */
                
                for (const auto &i : s2c)
                {
                    const auto m = r.match(i.first);
                    
                    const auto expected = m->concent();
                    auto measured = 0;
                    
                    auto x = 0.0;
                    auto y = 0.0;
                    
                    for (const auto &j : i.second)
                    {
                        x += c2m.at(j->id);
                        y += j->l.length();
                    }
                    
                    measured = x / y;
                    
                    stats.add(i.first, expected, measured);
                }
                
                //
                //                    assert(align.seq->l.length());
                //                    assert(contig.k_cov && contig.k_len);
                //
                //                    sumKLength += contig.k_len;
                //
                //                    // Normalized k-mer coverage
                //                    const auto n_cov = contig.normalized();
                //
                //                    switch (cov)
                //                    {
                //                        case WendySmooth:    { measured += n_cov * contig.k_len;          break; }
                //                        case KMerCov_Contig: { measured += n_cov / contig.k_len;          break; }
                //                        case KMerCov_Sequin: { measured += n_cov / align.seq->l.length(); break; }
                //                    }
                //
                //                    /*
                //                     * Calculate for the average depth for alignment and sequin
                //                     */
                //                    
                //                    align.depthAlign  += align.contigs[i].l.length() * contig.k_cov / align.contigs[i].l.length();
                //                    align.depthSequin += align.contigs[i].l.length() * contig.k_cov;
                
                
                
                break;
            }
        }
        
        /*
         * Post-processing... Add up the histogram...
         */
        
        switch (o.soft)
        {
            case MExpress::BWA:
            case MExpress::Bowtie:
            {
                for (auto &i : stats.hist)
                {
                    if (i.second)
                    {
                        stats.add(i.first, r.match(i.first)->concent(), i.second);
                    }
                    else
                    {
                        std::cout << i.first << std::endl;
                    }
                }
                
                break;
            }
                
            default : { break; }
        }
    }
    
    
/*
    switch (o.soft)
    {
        case Software::Velvet:  { stats.assembly = Velvet::analyze<MAssembly::Stats, DAsssembly::Contig>(file, &t);             break; }
        case Software::RayMeta: { stats.assembly = RayMeta::analyze<MAssembly::Stats, DAsssembly::Contig>(file, o.contigs, &t); break; }
    }
*/

    //stats.blat = t;
 
//    if (!stats.assembly.n)
//    {
//        throw std::runtime_error("No contig detected in the input file. Please check and try again.");
//    }
//    else if (stats.assembly.contigs.empty())
//    {
//        throw std::runtime_error("No contig aligned in the input file. Please check and try again.");
//    }

    //stats.n_chrT = stats.assembly.contigs.size();
    //stats.n_geno = stats.assembly.n - stats.n_chrT;

//    o.info("Analyzing the alignments");
//
//    for (auto &meta : stats.blat.metas)
//    {
//        auto &align = meta.second;
//        
//        if (!r.match(align->seq->id))
//        {
//            o.warn((boost::format("%1% not defined in the mixture. Skipped.") % align->seq->id).str());
//            continue;
//        }
//        
//        /*
//         * Calculate the limit of sensitivity. LOS is defined as the sequin with the lowest amount of
//         * concentration while still detectable in the experiment.
//         */
//
//        if (stats.limit.id.empty() || align->seq->concent() < stats.limit.abund)
//        {
//            stats.limit.id     = align->seq->id;
//            stats.limit.abund  = align->seq->concent();
//            stats.limit.counts = align->contigs.size();
//        }
//        
//        const auto p = MExpress::calculate(stats, stats.blat, stats.assembly, align->seq->id, *meta.second, o, o.coverage);
//
//        if (p.x && p.y)
//        {
//            stats.add(align->seq->id, p.x, p.y);
//        }
//    }

    stats.limit = r.absolute(stats.hist);

    return stats;
}

static void generateContigs(const FileName &file, const MExpress::Stats &stats, const MExpress::Options &o)
{
    o.info("Generating " + file);
    o.writer->open(file);

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.writer->write((boost::format(format) % "contigID"
                                           % "seqID"
                                           % "length"
                                           % "coverage"
                                           % "normalized").str());
    
    for (const auto &i : stats.blat.aligns)
    {
        if (stats.assembly.contigs.count(i.first))
        {
            const auto &contig = stats.assembly.contigs.at(i.first);
            
            o.writer->write((boost::format(format) % i.first
                                                   % i.second->id()
                                                   % contig.k_len
                                                   % contig.k_cov
                                                   % contig.normalized()).str());
        }
        else
        {
            o.writer->write((boost::format(format) % i.first
                                                   % i.second->id()
                                                   % "-"
                                                   % "-"
                                                   % "-").str());
        }
    }
    
    o.writer->close();
}

void MExpress::report(const std::vector<FileName> &files, const MExpress::Options &o)
{
    const auto stats = MExpress::analyze(files, o);

    /*
     * Generating MetaExpress_summary.stats
     */
    
    o.info("Generating MetaExpress_summary.stats");
    o.writer->open("MetaExpress_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rGeno,
                                                files[0],
                                                stats.hist,
                                                stats,
                                                stats,
                                                "sequins"));
    o.writer->close();
    
    /*
     * Generating MetaExpress_quins.stats
     */
    
    o.info("Generating MetaExpress_quins.stats");
    o.writer->open("MetaExpress_quins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating MetaExpress_express.R
     */
    
    o.info("Generating MetaExpress_express.R");
    o.writer->open("MetaExpress_express.R");
    o.writer->write(RWriter::createScript("MetaExpress_quins.stats", PlotMExpress()));
    o.writer->close();
    
    /*
     * Generating MetaExpress_contigs.stats
     */

    //generateContigs("MetaExpress_contigs.stats", stats, o);
}