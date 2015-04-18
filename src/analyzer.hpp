#ifndef GI_ANALYZER_HPP
#define GI_ANALYZER_HPP

#include <map>
#include <memory>
#include "types.hpp"
#include "writers/mock_writer.hpp"

namespace Spike
{
    #define INIT_COUNTER(x) std::map<SequinID, Counts> x; std::for_each(r.seqs_iA.begin(), r.seqs_iA.end(), [&](const std::pair<GeneID, Sequin> &p) { x[p.first] = 0; }); assert(x.size() == r.seqs_iA.size());

    #define ANALYZE_COUNTS(x,y) const auto y = Expression::analyze(x); assert(!y.limit_count || r.seqs_iA.count(y.limit_key));

    struct AnalyzerStats
    {
        // Total number of reads aligned
        Counts n = 0;
        
        // Number of reads aligned to the chromosome
        Counts nr = 0;
        
        // Number of reads aligned to the real sample
        Counts nq = 0;

        inline Percentage pr() const
        {
            return nr / n;
        }

        inline Percentage pq() const
        {
            return nq / n;
        }

        inline Percentage dilution() const
        {
            return nq ? static_cast<Percentage>(nr / nq) : 1;
        }
    };

    enum RNALevel
    {
        LevelGene,
        LevelIsoform,
    };

    template <typename Level> struct AnalyzerOptions
    {
        Level level;

        // How the results are written
        std::shared_ptr<Writer> writer = std::shared_ptr<Writer>(new MockWriter());
    };
}

#endif