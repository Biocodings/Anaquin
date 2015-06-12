#ifndef GI_ANALYZER_HPP
#define GI_ANALYZER_HPP

#include <map>
#include <memory>
#include "classify.hpp"
#include "data/types.hpp"
#include <boost/format.hpp>
#include <ss/regression/lm.hpp>
#include "writers/r_writer.hpp"
#include "stats/sensitivity.hpp"
#include "writers/mock_writer.hpp"

namespace Spike
{
    typedef std::map<Locus, Counts>     LocusCounter;
    typedef std::map<GeneID, Counts>    GeneCounter;
    typedef std::map<SequinID, Counts>  SequinCounter;
    typedef std::map<IsoformID, Counts> IsoformCounter;

    // Tracking for each sequin
    typedef std::map<SequinID, std::vector<Locus>> SequinTracker;

    template <typename Iter, typename T> static std::map<T, Counts> counter(const Iter &iter)
    {
        std::map<T, Counts> c;
        
        for (const auto &i : iter)
        {
            c[static_cast<T>(i)] = 0;
        }
        
        return c;
    }

    template <typename T> static void sums(const std::map<T, Counts> &m, Counts &c)
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

    struct Analyzer
    {
        template <typename T1, typename T2, typename Iter> static std::map<T1, T2> tracker(const Iter &iter)
        {
            std::map<T1, T2> m;

            for (const auto &i : iter)
            {
                m[static_cast<T1>(i)] = T2();
            }

            return m;
        }
    };

    struct DAnalyzer
    {
        static SequinCounter counterSequins()
        {
            return counter<std::vector<BedFeature>, SequinID>(Standard::instance().d_annot);
        }
    };

    class RAnalyzer
    {
        public:
            static SequinTracker sequinTracker()
            {
                return Analyzer::tracker<SequinID, std::vector<Locus>>(Standard::instance().r_sequins);
            }

            static LocusCounter exonCounter()
            {
                return counter<std::vector<Feature>, Locus>(Standard::instance().r_exons);
            }

            static LocusCounter intronCounter()
            {
                return counter<std::vector<Feature>, Locus>(Standard::instance().r_introns);
            }
        
            static GeneCounter geneCounter()
            {
                return counter<std::vector<Feature>, GeneID>(Standard::instance().r_genes);
            }
        
            static IsoformCounter isoformCounter()
            {
                return counter<std::vector<Sequin>, IsoformID>(Standard::instance().r_sequins);
            }
    };

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

    // Classify at the base-level by counting for non-overlapping regions
    template <typename I1, typename I2> void countBase(const I1 &r, const I2 &q, Confusion &m, Counter &c)
    {
        typedef typename I2::value_type Type;
        
        /*
         * The features in the query might overlap. We'd need to merge them before analysis.
         */
        
        const auto merged = Locus::merge<Type, Locus>(q);
        assert(!Locus::overlap(merged));

        for (const auto &l : merged)
        {
            m.nq   += l.length();
            m.tp() += countOverlaps(r, l, c);
            m.fp()  = m.nq - m.tp();
        }
    }

    struct ModelStats
    {
        Sensitivity s;

        std::vector<SequinID> z;
        
        // Known concentration for sequins
        std::vector<Concentration> x;
        
        // Measured coverage for sequins
        std::vector<Coverage> y;
        
        LinearModel linear() const
        {
            const auto m = SS::lm("y ~ x", SS::data.frame(SS::c(y), SS::c(x)));
            
            LinearModel lm;
            
            // Pearson correlation
            lm.r = SS::cor(x, y);
            
            // Adjusted R2
            lm.r2 = m.ar2;
            
            // Constant coefficient
            lm.c = m.coeffs[0].v;
            
            // Regression slope
            lm.m = m.coeffs[1].v;
            
            return lm;
        }
    };

    struct AnalyzerOptions
    {
        std::set<SequinID> filters;

        // By default, the results are written to no-where
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
        template <typename T> static void script(const std::string &file,
                                                 const std::vector<T> &x,
                                                 const std::vector<T> &y,
                                                 const std::vector<std::string> &z,
                                                 std::shared_ptr<Writer> writer)
        {
            writer->open(file);
            writer->write(RWriter::write(x, y, z));
            writer->close();
        }

        static void script(const std::string &file,
                           const ModelStats  &stats,
                           std::shared_ptr<Writer> writer)
        {
            writer->open(file);
            writer->write(RWriter::write(stats.x, stats.y, stats.z));
            writer->close();
        }

        template <typename ID, typename Stats, typename Writer>
            static void report(const std::string &name,
                               const std::string &r,
                               const Stats &stats,
                               const std::map<ID, Counts> &c,
                               std::shared_ptr<Writer> writer)
        {
            const std::string format = "%1%\t%2%\t%3%\t%4%";
            const auto &lm = stats.linear();
            
            writer->open(name);
            writer->write((boost::format(format) % "r" % "slope" % "r2" % "ss").str());
            writer->write((boost::format(format) % lm.r % lm.m % lm.r2 % stats.s.abund).str());
            writer->write("\n");

            for (const auto &p : c)
            {
                writer->write((boost::format("%1%\t%2%") % p.first % p.second).str());
            }
            
            writer->close();

            /*
             * Generate a plot for the fold-change relationship
             */

            writer->open(r);
            writer->write(RWriter::write(stats.x, stats.y, stats.z));
            writer->close();
        }

        template <typename Writer> static void report(const FileName &name,
                                                      const Performance &p,
                                                      const Counter &c,
                                                      Writer writer)
        {
            const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";

            writer->open(name);
            writer->write((boost::format(format) % "sn" % "sp" % "los" % "ss" % "counts").str());
            writer->write((boost::format(format) % p.m.sn()
                                                 % p.m.sp()
                                                 % p.s.id
                                                 % p.s.counts
                                                 % p.s.counts).str());
            writer->write("\n");

            for (const auto &p : c)
            {
                writer->write((boost::format("%1%\t%2%") % p.first % p.second).str());
            }

            writer->close();
        };
    };
}

#endif