#ifndef GI_ANALYZER_HPP
#define GI_ANALYZER_HPP

#include <map>
#include <memory>
#include "types.hpp"
#include "classify.hpp"
#include "sensitivity.hpp"
#include <boost/format.hpp>
#include "writers/mock_writer.hpp"

namespace Spike
{
    typedef std::map<std::string, Counts> Counter;

    struct LinearModel
    {
        // Least-squared constant coefficient
        double c;

        // Least-squared slope coefficient
        double m;
        
        // Adjusted R2
        double r2;

        // Pearson correlation
        double r;
    };
    
    template <typename Iter1, typename Iter2> void countBase(const Iter1 &r, const Iter2 &q, Confusion &m)
    {
        /*
         * Classify by constructing non-overlapping region for the query
         */

        const auto q_merged = Locus::merge(q);
        assert(!Locus::overlap(q_merged));

        for (auto l : q_merged)
        {
            m.nq   += l.length();
            m.tp() += countOverlaps(r, l);
            m.fp()  = m.nq - m.tp();
        }
    }
    
    inline std::map<SequinID, Counts> countsForSequins()
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

    inline std::map<GeneID, Counts> countsForGenes()
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
    
    class AnalyzerStats
    {
        public:
            //inline Percentage dilution() const
            //{
            //    return n() ? static_cast<Percentage>(nt) / n() : 1.0;
            //}
    };

    struct AnalyzerOptions
    {
        std::set<SequinID> filters;

        // How the results are written
        std::shared_ptr<Writer> writer = std::shared_ptr<Writer>(new MockWriter());
    };

    struct SingleMixtureOptions : public AnalyzerOptions
    {
        Mixture mix = MixA;
    };

    struct DoubleMixtureOptions : public AnalyzerOptions
    {
        const Mixture rMix = MixA;
        const Mixture qMix = MixB;
    };

    struct AnalyzeReporter
    {
        template <typename ID> static void report
                    (const std::string &name,
                     const Confusion &m,
                     const Sensitivity &s,
                     const std::map<ID, Counts> &c,
                     std::shared_ptr<Writer> writer)
        {
            const std::string format = "%1%\t%2%\t%3%";
            
            writer->open(name);
            writer->write((boost::format(format) % "sp" % "sn" % "ss").str());
            writer->write((boost::format(format) % m.sp() % m.sn() % s.abund).str());
            writer->write("\n");

            for (auto p : c)
            {
                writer->write((boost::format("%1%\t%2%") % p.first % p.second).str());
            }
            
            writer->close();
        }
    };
}

#endif