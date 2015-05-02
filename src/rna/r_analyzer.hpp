#ifndef GI_R_ANALYZER_HPP
#define GI_R_ANALYZER_HPP

#include "biology.hpp"
#include "analyzer.hpp"

namespace Spike
{
    struct RAnalyzer
    {
        protected:
            static Counter counter(RNALevel level, Mixture mix)
            {
                const auto &s = Standard::instance();
                std::map<SequinID, Counts> m;

                if (level == Gene)
                {
                    const auto p = s.r_pair(mix);

                    std::for_each(p.begin(), p.end(), [&](const std::pair<SequinID, Sequins> &p)
                                  {
                                      m[p.first] = 0;
                                  });
                    assert(m.size() == s.r_seqs_gA.size() && m.size() == s.r_seqs_gB.size());

                    return m;
                }
                else
                {
                    const auto se = s.r_sequin(mix);

                    std::for_each(se.begin(), se.end(), [&](const std::pair<SequinID, Sequin> &p)
                                  {
                                      m[p.first] = 0;
                                  });
                    assert(m.size() == s.r_seqs_iA.size() && m.size() == s.r_seqs_iB.size());

                    return m;
                }
            }
    };
}

#endif