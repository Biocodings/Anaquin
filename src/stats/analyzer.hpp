#ifndef GI_ANALYZER_HPP
#define GI_ANALYZER_HPP

#include <map>
#include <memory>
#include <boost/format.hpp>
#include "stats/classify.hpp"
#include <ss/regression/lm.hpp>
#include "writers/r_writer.hpp"
#include "stats/sensitivity.hpp"
#include "writers/mock_writer.hpp"

namespace Anaquin
{
    /*
     * List of softwares supported by Anaquin
     */

    enum Software
    {
        NA,
        Star,
        TopHat,
        Cufflink,
    };

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
        // Empty Implementation
    };

    struct MappingStats
    {
        // Total mapped to the in-silico chromosome
        Counts n_chrT = 0;

        // Total mapped to the human genome
        Counts n_hg38 = 0;

        // Fraction of sequin spiked
        inline Percentage dilution() const { return n_chrT / (n_chrT + n_hg38); }
    };
    
    /*
     * Represents a sequin that is not detected in the experiment
     */
    
    struct MissingSequin
    {
        MissingSequin(const SequinID &id, Concentration abund) : id(id), abund(abund) {}

        SequinID id;

        // The expect abundance
        Concentration abund;
    };

    typedef std::vector<MissingSequin> MissingSequins;

    /*
     * Represents a simple linear regression fitted by maximum-likehihood estimation.
     *
     *   Model: y ~ c + m*x
     */

    struct LinearModel
    {
        // Constant coefficient
        double c;

        // Least-squared slope coefficient
        double m;

        // Adjusted R2
        double r2;

        // Pearson correlation
        double r;
        
        // Adjusted R2
        double ar2;

        double f, p;
        double sst, ssm, sse;

        // Degree of freedoms
        unsigned sst_df, ssm_df, sse_df;
    };

    // Classify at the base-level by counting for non-overlapping regions
    template <typename I1, typename I2> void countBase(const I1 &r, const I2 &q, Confusion &m, SequinHist &c)
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

    struct Point
    {
        Point(double x = 0.0, double y = 0.0) : x(x), y(y) {}
        
        // Data point for the coordinate
        double x, y;
    };

    struct LinearStats : public std::map<SequinID, Point>
    {
        Sensitivity s;

        inline void add(const SequinID &id, double x, double y)
        {
            (*this)[id] = Point(x, y);
        }

        inline LinearModel linear() const
        {
            std::vector<double> x, y;
            
            for (const auto &p : *this)
            {
                if (!isnan(p.second.x) && !isnan(p.second.y))
                {
                    x.push_back(p.second.x);
                    y.push_back(p.second.y);
                }
            }

            const auto m = SS::lm("y~x", SS::data.frame(SS::c(y), SS::c(x)));

            LinearModel lm;
            
            lm.f   = m.f;
            lm.p   = m.p;
            lm.r2  = m.ar2;
            lm.ar2 = m.ar2;
            lm.r   = SS::cor(x, y);
            lm.sst = m.total.ss;
            lm.ssm = m.model.ss;
            lm.sse = m.error.ss;
            lm.c   = m.coeffs[0].value;
            lm.m   = m.coeffs[1].value;

            lm.sst_df = m.total.df;
            lm.ssm_df = m.model.df;
            lm.sse_df = m.error.df;

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

    struct FuzzyOptions : public AnalyzerOptions
    {
        double fuzzy;
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

    struct ViewerOptions : public AnalyzerOptions
    {
        std::string path;
    };

    struct Expression_
    {
        template <typename Ref> static Sensitivity calculate(const std::map<std::string, Counts> &h, const Ref &r)
        {
            Sensitivity s;
            
            // The lowest count must be zero because it can't be negative
            s.counts = std::numeric_limits<unsigned>::max();
            
            for (auto iter = h.begin(); iter != h.end(); iter++)
            {
                const auto counts = iter->second;
                
                /*
                 * Is this sequin detectable? If it's detectable, what about the concentration?
                 * By definition, detection limit is defined as the smallest abundance while
                 * still being detected.
                 */
                
                if (counts)
                {
                    const auto &id = iter->first;
                    const auto seq = r.findGene(id);
                    
                    // Hard to believe a sequin in the histogram is undefined
                    assert(seq);
                    
                    if (counts < s.counts || (counts == s.counts && seq->abund(MixA) < s.abund))
                    {
                        s.id     = id;
                        s.counts = counts;
                        s.abund  = seq->abund(MixA);
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
        template <typename Stats, typename Writer> static void missing(const std::string &file,
                                                                       const Stats &stats,
                                                                       Writer writer)
        {
            const auto format = "%1%\t%2%";

            writer->open(file);
            writer->write((boost::format(format) % "id" % "abund").str());

            for (const auto &i : stats.miss)
            {
                writer->write((boost::format(format) % i.id % i.abund).str());
            }

            writer->close();
        }

        template <typename Writer> static void writeCSV(const std::vector<double> &x,
                                                        const std::vector<double> &y,
                                                        const std::vector<std::string> &z,
                                                        const std::string &file,
                                                        Writer writer)
        {
            writer->open(file);
            writer->write("ID,expect,measure");

            /*
             * Prefer to write results in sorted order
             */

            std::set<std::string> sorted(z.begin(), z.end());

            for (const auto &s : sorted)
            {
                const auto it = std::find(z.begin(), z.end(), s);
                const auto i  = std::distance(z.begin(), it);

                writer->write((boost::format("%1%,%2%,%3%") % z.at(i) % x.at(i) % y.at(i)).str());
            }

            writer->close();
        }

        /*
         * Provides a common framework for generating a R scatter plot
         */

        template <typename Stats, typename Writer> static void scatter(const Stats &stats,
                                                                       const std::string prefix,
                                                                       const std::string unit,
                                                                       Writer writer)
        {
            //assert(stats.x.size() == stats.y.size() && stats.y.size() == stats.z.size());
            
            std::vector<double> x, y;
            std::vector<std::string> z;
            
            /*
             * Ignore any invalid value...
             */
            
            for (const auto &p : stats)
            {
                if (!isnan(p.second.x) && !isnan(p.second.y))
                {
                    z.push_back(p.first);
                    x.push_back(p.second.x);
                    y.push_back(p.second.y);
                }
            }
            
            /*
             * Generate a script for data visualization
             */
            
            writer->open(prefix + "_plot.R");
            writer->write(RWriter::write(x, y, z, unit, stats.s.abund));
            writer->close();
            
            /*
             * Generate CSV for each sequin
             */

            writeCSV(x, y, z, prefix + "_plot.csv", writer);
        }

        template <typename Stats, typename Writer> static void linear(const Stats &stats,
                                                                      const std::string prefix,
                                                                      const std::string unit,
                                                                      Writer writer,
                                                                      bool r       = true,
                                                                      bool summary = true,
                                                                      bool sequin  = true)
        {
            //assert(stats.x.size() == stats.y.size() && stats.y.size() == stats.z.size());

            std::vector<double> x, y;
            std::vector<std::string> z;

            /*
             * Ignore any invalid value...
             */

            for (const auto &p : stats)
            {
                if (!isnan(p.second.x) && !isnan(p.second.y))
                {
                    z.push_back(p.first);
                    x.push_back(p.second.x);
                    y.push_back(p.second.y);
                }
            }

            /*
             * Generate a script for data visualization
             */

            if (r)
            {
                writer->open(prefix + ".R");
                writer->write(RWriter::write(x, y, z, unit, stats.s.abund));
                writer->close();
            }
            
            const std::string format = "%1%\t%2%\t%3%\t%4%";
            const auto lm = stats.linear();

            /*
             * Generate summary statistics
             */

            if (summary)
            {
                writer->open(prefix + "_summary.stats");
                writer->write((boost::format(format) % "r" % "slope" % "r2" % "ss").str());
                writer->write((boost::format(format) % lm.r % lm.m % lm.r2 % stats.s.abund).str());
                writer->close();
            }

            /*
             * Generate CSV for each sequin
             */

            if (sequin)
            {
                //writeCSV(x, y , z, prefix + "_quins.csv", writer);
            }
        }
    };
}

#endif
