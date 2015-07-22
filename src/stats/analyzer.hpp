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

namespace Anaquin
{
    typedef std::map<Locus, Counts>    LocusCounter;
    typedef std::map<GeneID, Counts>   GeneCounter;
    typedef std::map<SequinID, Counts> SequinCounter;

    typedef std::map<Locus, Counts>    LocusHist;
    typedef std::map<SequinID, Counts> SequinHist;

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
        for (const auto &i : m)
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
        /*
         * Create a histogram with a key type of R. The key type is usually the SequinID, but can also
         * be GeneID or anything else. The histogram is initialized to zero.
         */

        template <typename T1, typename T2, typename R = SequinID> static std::map<R, Counts>
                    histogram(const std::map<T1, T2> &m)
        {
            std::map<R, Counts> c;

            for (const auto &p : m)
            {
                c[static_cast<R>(p.first)] = 0;
            }

            return c;
        }

        
        
        static std::map<Locus, Counts> histogram__(const std::map<Locus, Variation> &m)
        {
            std::map<Locus, Counts> c;
            
            for (const auto &p : m)
            {
                c[p.first] = 0;
            }
            
            return c;
        }
        
        
        
        
        
        
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

    class TAnalyzer
    {
        public:
            typedef std::map<SequinID, Performance> GeneTracker;

            static GeneTracker geneTracker()
            {
                return Analyzer::tracker<GeneID, Performance>(Standard::instance().r_genes);
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

            static SequinCounter sequinCounter()
            {
                return counter<std::set<SequinID>, SequinID>(Standard::instance().r_seqIDs);
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
        
        assert(!Locus::overlap(r));
        
        const auto merged = Locus::merge<Type, Locus>(q);
        
        for (const auto &l : merged)
        {
            m.nq   += l.length();
            m.tp() += countOverlaps(r, l, c);
            m.fp()  = m.nq - m.tp();
            
            // Make sure we don't run into negative
            assert(m.nq >= m.tp());
        }

        assert(!Locus::overlap(merged));
    }

    struct ModelStats
    {
        Sensitivity s;

        // Sequin IDs for each x and y
        std::vector<SequinID> z;

        // Known concentration for sequins
        std::vector<Concentration> x;

        // Measured coverage for sequins
        std::vector<Coverage> y;

        inline LinearModel linear() const
        {
            const auto m = SS::lm("y~x", SS::data.frame(SS::c(y), SS::c(x)));
            
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

        std::shared_ptr<Writer> writer = std::shared_ptr<Writer>(new MockWriter());
        std::shared_ptr<Writer> logger = std::shared_ptr<Writer>(new MockWriter());
        std::shared_ptr<Writer> output = std::shared_ptr<Writer>(new MockWriter());

        enum LogLevel
        {
            Info,
            Warn,
            Error,
        };

        inline void warn(const std::string &s) const
        {
            logger->write("[WARN]: " + s);
            output->write("[WARN]: " + s);
        }

        inline void wait(const std::string &s) const
        {
            logger->write("[WAIT]: " + s);
            output->write("[WAIT]: " + s);
        }
        
        inline void info(const std::string &s) const
        {
            logInfo(s);
            output->write("[INFO]: " + s);
        }

        inline void logInfo(const std::string &s) const
        {
            logger->write("[INFO]: " + s);
        }

        inline void logWarn(const std::string &s) const
        {
            logger->write("[WARN]: " + s);
        }
        
        inline void error(const std::string &s) const
        {
            logger->write("[ERROR]: " + s);
            output->write("[ERROR]: " + s);
        }

        // Write to the standard terminal
        inline void out(const std::string &s) const { output->write(s); }
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

    struct Expression
    {
        template <typename Map> static void print(const Map &m)
        {
            for (auto iter = m.begin(); iter != m.end(); iter++)
            {
                std::cout << iter->first << "  " << iter->second << std::endl;
            }
        }
        
        template <typename T, typename ID, typename S> static Sensitivity
            analyze(const std::map<T, Counts> &c, const std::map<ID, S> &m)
        {
            Sensitivity s;
            
            // The lowest count must be zero because it can't be negative
            s.counts = std::numeric_limits<unsigned>::max();

            for (auto iter = c.begin(); iter != c.end(); iter++)
            {
                const auto counts = iter->second;
                
                /*
                 * Is this sequin detectable? If it's detectable, what about the concentration?
                 * By definition, the detection limit is defined as the smallest abundance while
                 * still being detected.
                 */
                
                if (counts)
                {
                    const auto &id = iter->first;
                    
                    if (!m.count(id))
                    {
                        //std::cout << "Warning: " << id << " not found in LOS" << std::endl;
                        continue;
                    }
                    
                    if (counts < s.counts || (counts == s.counts && m.at(id).abund() < s.abund))
                    {
                        s.id     = id;
                        s.counts = counts;
                        s.abund  = m.at(s.id).abund();
                    }
                }
            }
            
            if (s.counts == std::numeric_limits<unsigned>::max())
            {
                s.counts = 0;
            }
            
            return s;
        }
    };

    struct AnalyzeReporter
    {
        /*
         * Write the following for a linear regression model:
         *
         *    - Linear plot
         *    - Linear statistics
         *    - Linear CSV (data points)
         */
        
        template <typename Stats, typename Writer> static void linear(const Stats &stats,
                                                                      const std::string prefix,
                                                                      const std::string unit,
                                                                      Writer writer)
        {
            assert(stats.x.size() == stats.y.size() && stats.y.size() == stats.z.size());

            const std::string format = "%1%\t%2%\t%3%\t%4%";
            const auto lm = stats.linear();

            /*
             * Generate linear statistics
             */

            writer->open(prefix + ".stats");
            writer->write((boost::format(format) % "r" % "slope" % "r2" % "ss").str());
            writer->write((boost::format(format) % lm.r % lm.m % lm.r2 % stats.s.abund).str());
            writer->close();
            
            /*
             * Generate linear CSV
             */
            
            writer->open(prefix + ".csv");
            writer->write("ID\expect\tmeasure");
            
            for (std::size_t i = 0; i < stats.x.size(); i++)
            {
                writer->write((boost::format("%1%\t%2%\t%3%") % stats.z[i] % stats.x[i] % stats.y[i]).str());
            }

            writer->close();

            /*
             * Generate a linear plot
             */
            
            writer->open(prefix + ".R");
            writer->write(RWriter::write(stats.x, stats.y, stats.z, unit, stats.s.abund));
            writer->close();
        }

        template <typename Writer, typename Histogram> static void stats(const FileName &name,
                                                                         const Performance &p,
                                                                         const Histogram   &h,
                                                                         Writer writer)
        {
            const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";

            const auto sn = p.m.sn();
            const auto sp = p.m.sp();
            const auto ss = p.s.abund;

            assert(ss >= 0);
            assert(isnan(sn) || (sn >= 0 && sn <= 1.0));
            assert(isnan(sp) || (sp >= 0 && sp <= 1.0));

            writer->open(name);
            writer->write((boost::format(format) % "sn" % "sp" % "los"  % "ss" % "counts").str());
            writer->write((boost::format(format) %  sn  %  sp  % p.s.id %  ss  % p.s.counts).str());
            writer->write("\n");

            for (const auto &p : h)
            {
                writer->write((boost::format("%1%\t%2%") % p.first % p.second).str());
            }

            writer->close();
        };
    };
}

#endif