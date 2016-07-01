#include "data/standard.hpp"
#include "MetaQuin/m_kabund.hpp"
#include "parsers/parser_sam.hpp"
#include "parsers/parser_quast.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotMKAbund();

MKAbund::Stats MKAbund::analyze(const std::vector<FileName> &files, const MKAbund::Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    MKAbund::Stats stats;
    stats.hist = r.hist();

    switch (o.soft)
    {
        case Software::Kallsito:
        {
            ParserKallisto::parse(Reader(files[0]), [&](const ParserKallisto::Data &d, const ParserProgress &)
            {
                const auto m = r.match(d.id);
                
                if (m)
                {
                    // Expected abundance
                    const auto known = m->mixes.at(Mix_1);
                    
                    // Measured abundance
                    const auto measured = d.abund;
                    
                    if (measured)
                    {
                        stats.add(d.id, known, measured);
                        stats.n_syn++;
                        //stats.hist.at(d.id)++;
                    }
                }
            });
            
            break;
        }
            
        default:
        {
            // Eg: Contigs.fasta
            const auto contigs = files[0];
            
            // Eg: Contigs.tsv
            const auto abund = files[1];
            
            // Eg: align.psl
            const auto align = files[2];
            
            o.analyze(contigs);
            
            std::map<ContigID, Base> c2l;
            std::map<ContigID, SequinID> c2s;
            std::map<SequinID, std::vector<ContigID>> s2c;
            
            switch (o.aligner)
            {
                case MAligner::Blat:
                {
                    const auto x = MBlat::analyze(files[1]);
                    
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
                    
                    break;
                }
                    
                case MAligner::MetaQuast:
                {
                    ParserQuast::parseAlign(Reader(align), [&](const ParserQuast::ContigData &x,
                                                               const ParserProgress &)
                    {
                        if (r.match(x.id))
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
                        }
                    });
                    
                    break;
                }
                    
                case MAligner::Kallisto: { break; }
            }
            
            assert(!c2s.empty());
            assert(!c2l.empty());
            assert(!s2c.empty());
            
            // Mapping from contigs to k-mer coverage
            std::map<ContigID, Coverage> c2m;
            
            // Mapping from contigs to k-mer length
            std::map<ContigID, Base> c2kl;
            
            switch (o.soft)
            {
                case Software::Kallsito:
                {
                    break;
                }
                    
                case Software::Velvet:
                {
                    break;
                }
                    
                case Software::RayMeta:
                {
                    ParserTSV::parse(Reader(abund), [&](const ParserTSV::TSV &x, const ParserProgress &)
                    {
                        c2m[x.id]  = x.kmer;
                        c2kl[x.id] = x.klen;
                    });
                    
                    assert(!c2m.empty());
                    
                    /*
                     * Quantifying k-mer abundance
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

            break;
        }
    }

//    stats.limit = r.detectLimit(stats.hist);

    return stats;
}

//static void generateContigs(const FileName &file, const MKAbund::Stats &stats, const MKAbund::Options &o)
//{
//    o.info("Generating " + file);
//    o.writer->open(file);
//
//    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
//    
//    o.writer->write((boost::format(format) % "contigID"
//                                           % "seqID"
//                                           % "length"
//                                           % "coverage"
//                                           % "normalized").str());
//    
//    for (const auto &i : stats.blat.aligns)
//    {
//        if (stats.assembly.contigs.count(i.first))
//        {
//            const auto &contig = stats.assembly.contigs.at(i.first);
//            
//            o.writer->write((boost::format(format) % i.first
//                                                   % i.second->id()
//                                                   % contig.k_len
//                                                   % contig.k_cov
//                                                   % contig.normalized()).str());
//        }
//        else
//        {
//            o.writer->write((boost::format(format) % i.first
//                                                   % i.second->id()
//                                                   % "-"
//                                                   % "-"
//                                                   % "-").str());
//        }
//    }
//    
//    o.writer->close();
//}

void MKAbund::report(const std::vector<FileName> &files, const MKAbund::Options &o)
{
    const auto stats = MKAbund::analyze(files, o);

    /*
     * Generating MetaKAbund_summary.stats
     */
    
    o.generate("MetaKAbund_summary.stats");
    o.writer->open("MetaKAbund_summary.stats");
//    o.writer->write(StatsWriter::inflectSummary(o.rAnnot,
//                                                o.rAnnot,
//                                                files[0],
//                                                stats.hist,
//                                                stats,
//                                                stats,
//                                                "sequins"));
    o.writer->close();
    
    /*
     * Generating MetaKAbund_quins.stats
     */
    
    o.generate("MetaKAbund_sequins.csv");
    o.writer->open("MetaKAbund_sequins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating MetaKAbund_kCounts.R
     */
    
    o.generate("MetaKAbund_kCounts.R");
    o.writer->open("MetaKAbund_kCounts.R");
    o.writer->write(RWriter::createScript("MetaKAbund_sequins.csv", PlotMKAbund()));
    o.writer->close();

    /*
     * Generating MetaKAbund_contigs.stats
     */

    //generateContigs("MetaKAbund_contigs.stats", stats, o);
}