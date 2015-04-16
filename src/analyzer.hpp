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
        // Percentage of reads aligned with the reference chromosome
        Percentage pr;
        
        // Percentage of reads aligned with the query samples
        Percentage pq;
        
        // Total number of reads aligned
        Reads n = 0;
        
        // Number of reads aligned to the chromosome
        Reads nr = 0;
        
        // Number of reads aligned to the real sample
        Reads nq = 0;
        
        inline Percentage dilution() const
        {
            return nq ? static_cast<Percentage>(nr / nq) : 1;
        }
    };

    template <typename Mode> struct AnalyzerOptions
    {
        Mode mode;

        // How the results are written
        std::shared_ptr<Writer> writer = std::shared_ptr<Writer>(new MockWriter());
    };
}

#endif