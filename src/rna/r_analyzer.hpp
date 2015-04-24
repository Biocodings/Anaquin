#ifndef GI_R_ANALYZER_HPP
#define GI_R_ANALYZER_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct RAnalyzer
    {
        protected:
        
            static std::map<SequinID, Counts> countsForSequins()
            {
                const auto &r = Standard::instance();
            
                std::map<TranscriptID, Counts> m;
                std::for_each(r.r_seqs_iA.begin(), r.r_seqs_iA.end(), [&](const std::pair<TranscriptID, Sequin> &p)
                {
                    m[p.first] = 0;
                });
            
                assert(m.size() != r.r_seqs_gA.size() && m.size() != r.r_seqs_gB.size());
                assert(m.size() == r.r_seqs_iA.size() && m.size() == r.r_seqs_iB.size());
            
                return m;
            }

            static std::map<GeneID, Counts> countsForGenes()
            {
                const auto &r = Standard::instance();
            
                std::map<GeneID, Counts> m;
                std::for_each(r.r_seqs_gA.begin(), r.r_seqs_gA.end(), [&](const std::pair<GeneID, Sequins> &p)
                {
                    m[p.first] = 0;
                });
            
                assert(m.size() == r.r_seqs_gA.size() && m.size() == r.r_seqs_gB.size());
                assert(m.size() != r.r_seqs_iA.size() && m.size() != r.r_seqs_iB.size());
            
                return m;
            }
    };
}

#endif