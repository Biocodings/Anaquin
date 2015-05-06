#ifndef GI_R_ANALYZER_HPP
#define GI_R_ANALYZER_HPP

#include "analyzer.hpp"

namespace Spike
{
    typedef std::map<Locus, Counts>     LocusCounter;
    typedef std::map<GeneID, Counts>    GeneCounter;
    typedef std::map<IsoformID, Counts> IsoformCounter;

    class RAnalyzer
    {
        public:

            static LocusCounter exonCounter()
            {
                return RAnalyzer::counter<std::vector<Feature>, Locus>(Standard::instance().r_exons);
            }

            static LocusCounter intronCounter()
            {
                return RAnalyzer::counter<std::vector<Feature>, Locus>(Standard::instance().r_introns);
            }
        
            static GeneCounter geneCounter()
            {
                return RAnalyzer::counter<std::vector<Feature>, GeneID>(Standard::instance().r_genes);
            }

            static IsoformCounter isoformCounter()
            {
                return RAnalyzer::counter<std::vector<Sequin>, IsoformID>(Standard::instance().r_sequins);
            }

        private:
        
            template <typename Iter, typename T> static std::map<T, Counts> counter(const Iter &iter)
            {
                std::map<T, Counts> c;

                for (const auto &i : iter)
                {
                    c[static_cast<T>(i)] = 0;
                }

                return c;
            }

        protected:

            template <typename T> static void count_ref(const std::map<T, Counts> &m, Counts &c)
            {
                for (const auto & i : m)
                {
                    if (i.second == 0)
                    {
                        c++;
                    }
                    else
                    {
                        c += i.second;
                    }
                }

                assert(c);
            }
    };
}

#endif