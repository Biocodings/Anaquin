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
    template <typename Stats> void adjustFN(const Stats &s, Confusion &m)
    {
        m.fn = s.nr - m.tp; // TODO: Revise this formula
    }

    #define FIX_FN(x, y) adjustFN(x, y)

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
    
    struct AnalyzerStats
    {
        // Binary classification at the base level
        Confusion mb;
        
        // Limit of sensitivity at the base level
        Sensitivity sb;

        // Total number of samples
        Counts n = 0;
        
        // Number of samples aligned to the chromosome
        Counts nr = 0;
        
        // Number of samples aligned to the real sample
        Counts nq = 0;

        inline Percentage pr() const { return nr / n; }
        inline Percentage pq() const { return nq / n; }

        inline Percentage dilution() const
        {
            return n ? static_cast<Percentage>(nr) / n : 1.0;
        }
    };

    template <typename Level> struct AnalyzerOptions
    {
        Level level;

        // How the results are written
        std::shared_ptr<Writer> writer = std::shared_ptr<Writer>(new MockWriter());
    };
    
    template <typename Level> struct SingleMixtureOptions : public AnalyzerOptions<Level>
    {
        Mixture mix = MixA;
    };

    struct AnalyzeReporter
    {
        template <typename ID> static void reportCounts(const std::string &name,
                                                        const std::map<ID, Counts> &m,
                                                        std::shared_ptr<Writer> writer)
        {
            const std::string format = "%1%\t%2%";
            
            writer->open(name);
            writer->write((boost::format(format) % "name" % "counts").str());
            
            for (auto p : m)
            {
                writer->write((boost::format(format) % p.first % p.second).str());
            }
            
            writer->close();
        }

        static void reportClassify
            (const std::string &name, Percentage dilution, const Confusion &m,
                const Sensitivity &s, std::shared_ptr<Writer> writer)
        {
            const std::string format = "%1%\t%2%\t%3%\t%4%";
            
            writer->open(name);
            writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
            writer->write((boost::format(format) % dilution
                                                 % m.sp()
                                                 % m.sn()
                                                 % s.abund).str());
            writer->close();
        }
    };
}

#endif