#include "parsers/parser_sam.hpp"
#include "MetaQuin/m_kexpress.hpp"
#include "parsers/parser_quast.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMExpress();

MKExpress::Stats MKExpress::analyze(const std::vector<FileName> &files, const MKExpress::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MKExpress::Stats stats;
    stats.hist = r.hist();
    
    const auto contigF = files[0];
    
    //for (auto &file : files)
    {
        o.analyze(contigF);
        
        switch (o.soft)
        {
            case MKExpress::Velvet:
            {
                break;
            }
                
            case MKExpress::RayMeta:
            {
                std::map<ContigID, Base> c2l;
                std::map<ContigID, SequinID> c2s;
                std::map<SequinID, std::vector<ContigID>> s2c;
                
                /*
                 * What alignment methods is being used?
                 */
                
                const auto align = files[2];
                const auto isPSL = align.find("psl") != std::string::npos;
                
                if (isPSL)
                {
                    const auto x = MBlat::analyze(align);
                    
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
                }
                else
                {
                    ParserQuast::parseAlign(Reader(files[2]), [&](const ParserQuast::ContigData &x,
                                                                  const ParserProgress &)
                    {
                        for (const auto &c : x.contigs)
                        {
                            // Contigs.fasta doesn't have "_"
                            auto t = c;
                            
                            // Eg: contig-1056000000 2818 nucleotides
                            boost::replace_all(t, "_", " ");
                            
                            std::vector<std::string> toks;
                            Tokens::split(t, " ", toks);
                            
                            c2s[toks[0]] = x.id;
                            c2l[toks[0]] = stod(toks[1]);
                            s2c[x.id].push_back(toks[0]);
                        }
                    });
                }
                
                assert(!c2s.empty());
                assert(!c2l.empty());
                assert(!s2c.empty());
                
                // Mapping from contigs to k-mer coverage
                std::map<ContigID, Coverage> c2m;
                
                // Mapping from contigs to k-mer length
                std::map<ContigID, Base> c2kl;

                ParserTSV::parse(Reader(files[1]), [&](const ParserTSV::TSV &x, const ParserProgress &)
                {
                    c2m[x.id]  = x.kmer;
                    c2kl[x.id] = x.klen;
                });
                
                assert(!c2m.empty());

                /*
                 * Quantifying expressions
                 */
                
                for (const auto &i : s2c)
                {
                    const auto m = r.match(i.first);
                    
                    const auto expected = m->concent();
                    auto measured = 0.0;
                    
                    auto x = 0.0;
                    auto y = 0.0;
                    
                    for (const auto &j : i.second)
                    {
                        // K-mer length
                        const auto kl = c2kl.at(j);
                        
                        // Entire size of the contig (including bases that are not mapped)
                        const auto l = c2l[j];
                        
                        /*
                         * How should we normalize the k-mer observations? Should we normalize by the k-mer length?
                         * Should we normalize by the size of the contig?
                         *
                         * TODO: Is the k-mer observations reported in RayMeta already normalized?
                         */
                        
                        x += (double)c2m.at(j); // / kl;// / l;
                        y += l; //j->l.length();
                    }
                    
                    //measured = x / y;
                    measured = x;
                    
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
                //                     * Calculate the average depth for alignment and sequin
                //                     */
                //                    
                //                    align.depthAlign  += align.contigs[i].l.length() * contig.k_cov / align.contigs[i].l.length();
                //                    align.depthSequin += align.contigs[i].l.length() * contig.k_cov;
                
                
                
                break;
            }
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

static void generateContigs(const FileName &file, const MKExpress::Stats &stats, const MKExpress::Options &o)
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

void MKExpress::report(const std::vector<FileName> &files, const MKExpress::Options &o)
{
    const auto stats = MKExpress::analyze(files, o);

    /*
     * Generating MetaExpress_summary.stats
     */
    
    o.generate("MetaExpress_summary.stats");
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
    
    o.generate("MetaExpress_quins.stats");
    o.writer->open("MetaExpress_quins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating MetaExpress_express.R
     */
    
    o.generate("MetaExpress_express.R");
    o.writer->open("MetaExpress_express.R");
    o.writer->write(RWriter::createScript("MetaExpress_quins.stats", PlotMExpress()));
    o.writer->close();
    
    /*
     * Generating MetaExpress_contigs.stats
     */

    //generateContigs("MetaExpress_contigs.stats", stats, o);
}