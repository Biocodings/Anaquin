#ifndef GI_M_ASSEMBLY_HPP
#define GI_M_ASSEMBLY_HPP

#include <numeric>
#include "data/tokens.hpp"
#include "meta/m_blast.hpp"
#include "meta/m_histogram.h"
#include "stats/analyzer.hpp"
#include "parsers/parser_fa.hpp"

namespace Anaquin
{
    struct Contig
    {
        ContigID id;
        
        // The sequence being assembled
        Sequence seq;
        
        // Coverage in k-mer
        Coverage k_cov;
    };
    
    struct DAsssembly
    {
        template <typename T> struct Stats
        {
            Base mean, min, max;
            Base N20, N50, N80;

            // Total number of bases in contigs
            Base total;
            
            // Total number of bases in the assembly
            Base sum;

            // List of assembled contigs
            std::map<ContigID, T> contigs;
        };

        template <typename C = Contig, typename T = DAsssembly::Stats<C>> static T parse(const std::string &file, std::function<void (C&)> f)
        {
            T stats;
            Histogram h;

            /*
             * Read a file of assembled contigs. The file is assumed to be in a FA format.
             */
            
            ParserFA::parse(file, [&](const FALine &l, const ParserProgress &)
            {
                C c;

                c.id = l.id;

                // Sequence of the config
                c.seq = l.seq;

                // The histogram needs the length of the sequence
                h.insert(l.seq.size());

                // Allows to apply operations on a specific assembler
                f(c);

                stats.contigs[c.id] = c;
            });

            // This is copied from printContiguityStats() in Histogram.h of the Abyss source code.
            h = h.trimLow(500);
            
            /*
             * https://github.com/bcgsc/abyss/blob/e58e5a6666e0de0e6bdc15c81fe488f5d83085d1/Common/Histogram.h
             */
            
            stats.sum   = h.sum();
            stats.N50   = h.n50();
            stats.min   = h.minimum();
            stats.max   = h.maximum();
            stats.mean  = h.expectedValue();
            stats.N80   = h.weightedPercentile(1 - 0.8);
            stats.N20   = h.weightedPercentile(1 - 0.2);
            stats.sum   = h.sum();
            stats.total = std::accumulate(stats.contigs.begin(), stats.contigs.end(), 0,
                            [&](int sum, const std::pair<std::string, C> &p)
                            {
                                return sum + p.second.seq.size();
                            });
            return stats;
        }
    };

    struct Velvet
    {
        template <typename Stats, typename C> static Stats parse(const std::string &file)
        {
            Stats stats;

            /*
             * Read coverage from the contig file. The format looks like:
             *
             *      >NODE_77460_length_31_cov_1.129032
             */
            
            std::vector<std::string> tokens;

            return DAsssembly::parse<C, Stats>(file, [&](C &node)
            {
                Tokens::split(node.id, "_", tokens);

                // Parse the k-mer coverage
                node.k_cov = stod(tokens[tokens.size() - 1]);
            });
        }
    };

    struct MAssembly
    {
        enum Assembler
        {
            Velvet,
        };

        struct Stats : public DAsssembly::Stats<Contig>
        {
            // Statistics for the blat alignment
            MBlast::Stats blat;

            // Statistics for abundance
            LinearStats lm;
        };

        struct Options : public SingleMixtureOptions
        {
            Options() {}

            // Alignment file by blat
            FileName psl;

            // The type of the assembler used
            Assembler tool = Assembler::Velvet;
        };

        static Stats analyze(const std::string &, const Options &options = Options());        
        static Stats report (const std::string &, const Options &options = Options());
    };
}

#endif