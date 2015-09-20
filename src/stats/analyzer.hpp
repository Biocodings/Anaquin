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

    inline std::size_t countHist(const std::map<std::string, Counts> &m)
    {
        return std::count_if(m.begin(), m.end(), [&](const std::pair<SequinID, Counts> &i)
                                                 {
                                                     return i.second;
                                                 });
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
        inline Percentage dilution() const { return static_cast<double>(n_chrT) / (n_chrT + n_hg38); }
    };

    struct AlignmentStats : public MappingStats
    {
        Counts unmapped = 0;
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

    struct FusionStats : public MappingStats
    {
        // Number of fusions spanning across the genome and the synthetic chromosome
        Counts hg38_chrT;
    };

    struct LinearStats : public std::map<SequinID, Point>
    {
        Sensitivity s;

        inline void add(const SequinID &id, double x, double y)
        {
            (*this)[id] = Point(x, y);
        }

        /*
         * Compute a simple linear regression model. By default, this function assumes
         * the values are raw and will attempt to log-transform.
         */

        inline LinearModel linear(bool shouldLog = true) const
        {
            std::vector<double> x, y;

            auto f = [&](double v)
            {
                return shouldLog ? (v ? log2(v) : 0) : v;
            };

            for (const auto &p : *this)
            {
                if (!isnan(p.second.x) && !isnan(p.second.y))
                {
                    x.push_back(f(p.second.x));
                    y.push_back(f(p.second.y));
                }
            }

            const auto m = SS::lm("y~x", SS::R::data.frame(SS::R::c(y), SS::R::c(x)));

            LinearModel lm;
            
            lm.f      = m.f;
            lm.p      = m.p;
            lm.r2     = m.r2;
            lm.ar2    = m.ar2;
            lm.r      = SS::cor(x, y);
            lm.sst    = m.total.ss;
            lm.ssm    = m.model.ss;
            lm.sse    = m.error.ss;
            lm.c      = m.coeffs[0].value;
            lm.m      = m.coeffs[1].value;
            lm.sst_df = m.total.df;
            lm.ssm_df = m.model.df;
            lm.sse_df = m.error.df;

            return lm;
        }
    };

    struct WriterOptions
    {
        enum LogLevel
        {
            Info,
            Warn,
            Error,
        };

        std::shared_ptr<Writer> writer = std::shared_ptr<Writer>(new MockWriter());
        std::shared_ptr<Writer> logger = std::shared_ptr<Writer>(new MockWriter());
        std::shared_ptr<Writer> output = std::shared_ptr<Writer>(new MockWriter());

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

    struct AnalyzerOptions : public WriterOptions
    {
        std::set<SequinID> filters;
    };

    struct FuzzyOptions : public AnalyzerOptions
    {
        double fuzzy;
    };

    struct ViewerOptions : public AnalyzerOptions
    {
        std::string path;
    };

    struct SingleMixtureOptions : public AnalyzerOptions
    {
        Mixture mix = Mix_1;
    };

    struct DoubleMixtureOptions : public AnalyzerOptions
    {
        Mixture mix_1 = Mix_1;
        Mixture mix_2 = Mix_2;
    };

    struct AnalyzeReporter
    {
        /*
         * Provides a common framework to report a simple linear regression model
         */
        
        template <typename Writer, typename Stats> static void linear
            (const FileName &file, const Stats &stats, const std::string &units, Writer writer)
        {
            const auto summary = "Summary for dataset: %1%\n\n"
                                 "   Genome:      %2% %19%\n"
                                 "   Query:       %3% %19%\n"
                                 "   Reference:   %4% %19%\n\n"
                                 "   Sensitivity: %5% (attomol/ul) (%6%)\n"
                                 "   Detected:    %7% %19%\n\n"
                                 "   ***\n"
                                 "   *** The following statistics are computed on the log2 scale\n"
                                 "   ***\n"
                                 "   ***   Eg: If the data points are (1,1), (2,2). The correlation will\n"
                                 "   ***       be computed on (log2(1), log2(1)), (log2(2), log2(2)))\n"
                                 "   ***\n\n"
                                 "   Correlation:\t%8%\n"
                                 "   Slope:\t%9%\n"
                                 "   R2:\t%10%\n"
                                 "   F-statistic:\t%11%\n"
                                 "   P-value:\t%12%\n"
                                 "   SSM: %13%, DF: %14%\n"
                                 "   SSE: %15%, DF: %16%\n"
                                 "   SST: %17%, DF: %18%\n";
            const auto lm = stats.linear();
            
            writer->open(file);
            writer->write((boost::format(summary) % file
                                                  % stats.n_hg38
                                                  % stats.n_chrT
                                                  % stats.h.size()
                                                  % stats.ss.abund
                                                  % stats.ss.id
                                                  % countHist(stats.h)
                                                  % lm.r
                                                  % lm.m
                                                  % lm.r2
                                                  % lm.f
                                                  % lm.p
                                                  % lm.ssm
                                                  % lm.ssm_df
                                                  % lm.sse
                                                  % lm.sse_df
                                                  % lm.sst
                                                  % lm.sst_df
                                                  % units).str());
            writer->close();
        }
        
        template <typename Stats, typename Writer> static void missing(const FileName &file,
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
                                                                       const std::string &prefix,
                                                                       const std::string &xLabel,
                                                                       const std::string &yLabel,
                                                                       Writer writer)
        {
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
            writer->write(RWriter::coverage(x, y, z, xLabel, yLabel, stats.s.abund));
            writer->close();
            
            /*
             * Generate CSV for each sequin
             */

            writeCSV(x, y, z, prefix + "_plot.csv", writer);
        }
    };
}

#endif
