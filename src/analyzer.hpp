#ifndef GI_ANALYZER_HPP
#define GI_ANALYZER_HPP

#include <set>
#include <map>
#include <memory>
#include "types.hpp"
#include "classify.hpp"
#include "sensitivity.hpp"
#include <boost/format.hpp>
#include "writers/mock_writer.hpp"

namespace Spike
{
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
            // Classification at the base level
            Confusion mb;

            // Limit of sensitivity at the base level
            Sensitivity sb;

            // Counts for the reference
            Counts nr = 0;

            // Counts for the experiment
            Counts nq = 0;

            inline Counts n() const { return nr + nq; }

            inline Percentage pr() const { return nr / n(); }
            inline Percentage pq() const { return nq / n(); }

            inline Percentage dilution() const
            {
                return n() ? static_cast<Percentage>(nr) / n() : 1.0;
            }
    };

    template <typename Level> struct AnalyzerOptions
    {
        Level level;

        std::set<SequinID> filters;

        // How the results are written
        std::shared_ptr<Writer> writer = std::shared_ptr<Writer>(new MockWriter());
    };
    
    template <typename Level> struct SingleMixtureOptions : public AnalyzerOptions<Level>
    {
        Mixture mix = MixA;
    };

    struct AnalyzeReporter
    {
        template <typename ID> static void reportClassify
                    (const std::string &name,
                     Percentage dilution,
                     const Confusion &m,
                     const Sensitivity &s,
                     const std::map<ID, Counts> &c,
                     std::shared_ptr<Writer> writer)
        {
            const std::string format = "%1%\t%2%\t%3%\t%4%\n\n";
            
            writer->open(name);
            writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
            writer->write((boost::format(format) % dilution
                                                 % m.sp()
                                                 % m.sn()
                                                 % s.abund).str());

            for (auto p : c)
            {
                writer->write((boost::format("%1%\t%2%") % p.first % p.second).str());
            }
            
            writer->close();
        }
    };
}

#endif