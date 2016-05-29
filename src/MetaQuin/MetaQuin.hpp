#ifndef META_QUIN_HPP
#define META_QUIN_HPP

#include "MetaQuin/m_blat.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_tsv.hpp"
#include "MetaQuin/m_histogram.h"

namespace Anaquin
{
    enum class MAssembler
    {
        Velvet,
        RayMeta,
        Kallsito
    };
    
    enum class MAligner
    {
        MetaQuast,
        Blat,
        Kallisto,
    };

    struct DAsssembly
    {
        template <typename T> struct Stats : public MappingStats
        {
            Base mean, min, max;
            Base N20, N50, N80;
            
            // Total number of bases in contigs
            Base total;
            
            // Total number of bases in the assembly
            Base sum;
            
            // List of aligned contigs (it's not all contigs)
            std::map<ContigID, T> contigs;

            // Total number of contigs, whether it's aligned or not
            Counts n = 0;
        };

        struct Contig
        {
            ContigID id;
            
            // Length of the sequence being assembled
            Base len;
            
            // Size of the contig in k-mer
            Base k_len;
            
            // Unnormalized k-mer coverage
            Coverage k_cov;
            
            // Normalized k-mer coverage
            inline Coverage normalized() const { return k_cov ; }
        };

        struct DenoAssemblyImpl
        {
            virtual SequinID findSeq(const ContigID &) const = 0;
        };

        template <typename T = DAsssembly::Stats<Contig>> static T analyze
                    (const FileName &file, DenoAssemblyImpl *impl)
        {
            T stats;
            Histogram hist;

            ParserFA::parse(file, [&](const ParserFA::Data &x, const ParserProgress &)
            {
                Contig c;
                
                // Can this contig be mapped to a sequin?
                const auto seqID = impl->findSeq(c.id = x.id);
                
                if (seqID.empty())
                {
                    stats.n_geno++;
                    
                    // Don't bother if the contig isn't part of MetaQuins...
                    return;
                }
                else
                {
                    stats.n_chrT++;
                }
                
                // Size of the config
                c.len = x.seq.size();
                
                // The histogram needs the size of the sequence
                hist.insert(x.seq.size());
                
                stats.contigs[seqID] = c;

                /*
                if (blat)
                {
                    // Eg: "contig-30000000	21269" to "contig-30000000"
                    const auto first = Tokens::first(id, " ");
                    
                    id = blat->aligns.count(first) ? first : "";
                    
                    // Don't bother if the contig isn't defined in the alignment...
                    if (id.empty())
                    {
                        return;
                    }
                }
                 */
            });
            
            /*
             * https://github.com/bcgsc/abyss/blob/e58e5a6666e0de0e6bdc15c81fe488f5d83085d1/Common/Histogram.h
             */
            
            stats.sum   = hist.sum();
            stats.N50   = hist.n50();
            stats.min   = hist.minimum();
            stats.max   = hist.maximum();
            stats.mean  = hist.expectedValue();
            stats.N80   = hist.weightedPercentile(1 - 0.8);
            stats.N20   = hist.weightedPercentile(1 - 0.2);
            stats.sum   = hist.sum();
            stats.total = std::accumulate(stats.contigs.begin(), stats.contigs.end(), 0,
                                          [&](int sum, const std::pair<std::string, Contig> &p)
                                          {
                                              return sum + p.second.len;
                                          });
            return stats;
        }

        /*
         * Analyze a de-novo assembly given the assembly and blat alignment.
         */

        template <typename C = Contig, typename T = DAsssembly::Stats<C>> static T analyze_
                (const FileName &file, const MBlat::Stats *blat, std::function<void (C&)> f)
        {
            T stats;
            Histogram hist;

            ParserFA::parse(file, [&](const ParserFA::Data &data, const ParserProgress &)
            {
                //stats.n++;
                
                C c;
                c.id = data.id;
                
                auto id = c.id;
                
                if (blat)
                {
                    // Eg: "contig-30000000	21269" to "contig-30000000"
                    const auto first = Tokens::first(id, " ");

                    id = blat->aligns.count(first) ? first : "";

                    // Don't bother if the contig isn't defined in the alignment...
                    if (id.empty())
                    {
                        return;
                    }
                }
                
                // Size of the config
                c.len = data.seq.size();
                
                // The histogram needs the size of the sequence
                hist.insert(data.seq.size());
                
                // Allows to apply custom operation
                f(c);
                
                //stats.contigs[id] = c;
            });
            
            /*
             * https://github.com/bcgsc/abyss/blob/e58e5a6666e0de0e6bdc15c81fe488f5d83085d1/Common/Histogram.h
             */
            
//            stats.sum   = hist.sum();
//            stats.N50   = hist.n50();
//            stats.min   = hist.minimum();
//            stats.max   = hist.maximum();
//            stats.mean  = hist.expectedValue();
//            stats.N80   = hist.weightedPercentile(1 - 0.8);
//            stats.N20   = hist.weightedPercentile(1 - 0.2);
//            stats.sum   = hist.sum();
//            stats.total = std::accumulate(stats.contigs.begin(), stats.contigs.end(), 0,
//                            [&](int sum, const std::pair<std::string, C> &p)
//                            {
//                                return sum + p.second.len;
//                            });
            return stats;
        }
    };

    struct Velvet
    {
        template <typename Stats, typename C> static Stats analyze(const FileName &file,
                                                                   const MBlat::Stats *align = NULL)
        {
            Stats stats;
            
            /*
             * Read coverage from the contig file. The format looks like:
             *
             *      >NODE_77460_length_31_cov_1.129032
             */
            
            std::vector<std::string> toks;

            return DAsssembly::analyze_<C, Stats>(file, align, [&](C &node)
            {
                Tokens::split(node.id, "_", toks);

                node.k_len = stoi(toks[3]);
                node.k_cov = stod(toks[toks.size() - 1]) * node.k_len;
            });
        }
    };
    
    struct RayMeta
    {
        template <typename Stats, typename C> static Stats analyze(const FileName &file,
                                                                   const FileName &contigs,
                                                                   const MBlat::Stats *align)
        {
            std::map<ContigID, KMers> covs, lens;
            
            if (!contigs.empty())
            {
                ParserTSV::parse(Reader(contigs), [&](const ParserTSV::TSV &t, const ParserProgress &)
                {
                    covs[t.id] = t.kmer;
                    lens[t.id] = t.klen;
                });
            }

            Stats stats;
            
            return DAsssembly::analyze<C, Stats>(file, align, [&](C &node)
            {
                const auto id = Tokens::first(node.id, " ");
                
                if (covs.count(id))
                {
                    node.k_len = lens.at(id);
                    node.k_cov = covs.at(id);
                }
            });
        }
    };
}

#endif