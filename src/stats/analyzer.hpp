#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <map>
#include <memory>
#include "stats/limit.hpp"
#include <boost/format.hpp>
#include "stats/classify.hpp"
#include "writers/r_writer.hpp"
#include "writers/mock_writer.hpp"
#include <ss/regression/linear.hpp>

namespace Anaquin
{
    /*
     * Histogram manipulations and operations
     */

    template <typename T> static Counts sums(const std::map<T, Counts> &m)
    {
        Counts c = 0;

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

    // Number of elements in the histogram with at least an entry
    inline Counts detectHist(const std::map<std::string, Counts> &m)
    {
        return std::count_if(m.begin(), m.end(), [&](const std::pair<std::string, Counts> &i)
        {
            return i.second ? 1 : 0;
        });
    }
    
    struct Analyzer
    {
        // Empty Implementation
    };
    
    struct SingleInputStats
    {
        FileName src;
    };

    struct MappingStats
    {
        // The distribution of counts across chromosomes
        std::map<ChromoID, Counts> hist;

        // Total mapped to the in-silico chromosome
        Counts n_chrT = 0;

        // Total mapped to the experiment
        Counts n_expT = 0;

        // Fraction of sequin spiked
        inline Percentage dilution() const
        {
            return (n_chrT + n_expT) ? static_cast<double>(n_chrT) / (n_chrT + n_expT) : NAN;
        }
    };

    struct AlignmentStats : public MappingStats, public SingleInputStats
    {
        Counts unmapped = 0;

        template <typename T, typename F> void update(const T &t, F f)
        {
            if (!t.i)
            {
                if      (!t.mapped) { unmapped++; }
                else if (!f(t))     { n_expT++;   }
                else                { n_chrT++;   }
            }
        }
        
        template <typename T> void update(const T &t)
        {
            if (!t.i)
            {
                if      (!t.mapped)                       { unmapped++; }
                else if (t.id != Standard::instance().id) { n_expT++;   }
                else                                      { n_chrT++;   }
            }
        }
    };

//    /*
//     * Represents a sequin that is not detected in the experiment
//     */
//    
//    struct MissingSequin
//    {
//        MissingSequin(const SequinID &id, Concentration abund) : id(id), abund(abund) {}
//
//        SequinID id;
//
//        // The expect abundance
//        Concentration abund;
//    };
//
//    typedef std::vector<MissingSequin> MissingSequins;

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
        Counts hg38_chrT = 0;
    };

    struct LinearStats : public std::map<SequinID, Point>
    {
        Limit s;

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
                // Don't log zero...
                assert(!shouldLog || v);

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

            const auto m = SS::lm(SS::R::data.frame(SS::R::c(y), SS::R::c(x)));

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

        // Working directory
        FilePath working;
        
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
        
        inline void analyze(const std::string &s) const
        {
            info("Analyzing: " + s);
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

    struct SingleMixtureOption : public AnalyzerOptions
    {
        Mixture mix = Mix_1;
    };

    struct FuzzyOptions : public AnalyzerOptions
    {
        double fuzzy;
    };

    struct ViewerOptions : public AnalyzerOptions
    {
        Path path;
    };

    struct DoubleMixtureOptions : public AnalyzerOptions
    {
        Mixture mix_1 = Mix_1;
        Mixture mix_2 = Mix_2;
    };

    struct AnalyzeReporter
    {
        /*
         * Provides a common framework to report a linear regression model for two samples, typically
         * mixture A and B.
         */

        template <typename Writer, typename Stats_1, typename Stats_2, typename Stats> static void linear
                        (const FileName &f,
                         const FileName &d1,
                         const FileName &d2,
                         const Stats_1  &s1,
                         const Stats_2  &s2,
                         const Stats    &s,
                         const Units    &units,
                         Writer writer,
                         const Units &ref = "",
                         const Label &samples = "")
        {
            const auto summary = "Summary for dataset: %1% and %2%\n\n"
                                 "   %3% (A):       %4% %5%\n"
                                 "   Query (A):     %6% %5%\n\n"
                                 "   %3% (B):       %7% %5%\n"
                                 "   Query (B):     %8% %5%\n\n"
                                 "   Reference (A):   %9% %10%\n"
                                 "   Reference (B):   %11% %10%\n\n"
                                 "   Detected A:    %12% %10%\n"
                                 "   Detected B:    %13% %10%\n\n"
                                 "   Limit A:    %14% %15%\n"
                                 "   Limit B:    %16% %17%\n\n"
                                 "   Correlation: %18%\n"
                                 "   Slope:       %19%\n"
                                 "   R2:          %20%\n"
                                 "   F-statistic: %21%\n"
                                 "   P-value:     %22%\n"
                                 "   SSM:         %23%, DF: %24%\n"
                                 "   SSE:         %25%, DF: %26%\n"
                                 "   SST:         %27%, DF: %28%\n\n"
                                 "   ***\n"
                                 "   *** The following statistics are computed on the log2 scale.\n"
                                 "   ***\n"
                                 "   ***   Eg: If the data points are (1,1), (2,2). The correlation will\n"
                                 "   ***       be computed on (log2(1), log2(1)), (log2(2), log2(2)))\n"
                                 "   ***\n\n"
                                 "   Correlation: %29%\n"
                                 "   Slope:       %30%\n"
                                 "   R2:          %31%\n"
                                 "   F-statistic: %32%\n"
                                 "   P-value:     %33%\n"
                                 "   SSM:         %34%, DF: %35%\n"
                                 "   SSE:         %36%, DF: %37%\n"
                                 "   SST:         %38%, DF: %39%\n";

            const auto n_lm = s.linear(false);
            const auto l_lm = s.linear(true);
            
            writer->open(f);
            writer->write((boost::format(summary) % d1                          // 1
                                                  % d2
                                                  % (samples.empty() ? "Genome" : samples)
                                                  % s1.n_expT
                                                  % units
                                                  % s1.n_chrT
                                                  % s2.n_expT
                                                  % s2.n_chrT
                                                  % s1.h.size()
                                                  % (ref.empty() ? units : ref) // 10
                                                  % s2.h.size()
                                                  % detectHist(s1.h)            // 12
                                                  % detectHist(s2.h)            // 13
                                                  % s1.ss.abund                 // 14
                                                  % s1.ss.id                    // 15
                                                  % s2.ss.abund                 // 16
                                                  % s2.ss.id                    // 17
                                                  % n_lm.r                      // 18
                                                  % n_lm.m
                                                  % n_lm.r2
                                                  % n_lm.f
                                                  % n_lm.p                      // 22
                                                  % n_lm.ssm                    // 23
                                                  % n_lm.ssm_df
                                                  % n_lm.sse
                                                  % n_lm.sse_df
                                                  % n_lm.sst
                                                  % n_lm.sst_df                 // 28
                                                  % l_lm.r                      // 29
                                                  % l_lm.m
                                                  % l_lm.r2
                                                  % l_lm.f
                                                  % l_lm.p                      // 32
                                                  % l_lm.ssm                    // 33
                                                  % l_lm.ssm_df
                                                  % l_lm.sse
                                                  % l_lm.sse_df
                                                  % l_lm.sst
                                                  % l_lm.sst_df                 // 38
                           ).str());
            writer->close();
        }

        /*
         * Provides a common framework to report a simple linear regression model
         */

        template <typename Writer, typename Stats> static void linear
                        (const FileName &file,
                         const FileName &data,
                         const Stats &stats,
                         const Units &units,
                         Writer writer,
                         const Units &ref = "",
                         const Label &samples = "")
        {
            const auto summary = "Summary for dataset: %1%\n\n"
                                 "   %3%:         %2% %5%\n"
                                 "   Query:       %4% %5%\n\n"
                                 "   Reference:   %6% %7%\n"
                                 "   Detected:    %10% %7%\n\n"
                                 "   Sensitivity: %8% (attomol/ul) (%9%)\n\n"
                                 "   Correlation: %11%\n"
                                 "   Slope:       %12%\n"
                                 "   R2:          %13%\n"
                                 "   F-statistic: %14%\n"
                                 "   P-value:     %15%\n"
                                 "   SSM:         %16%, DF: %17%\n"
                                 "   SSE:         %18%, DF: %19%\n"
                                 "   SST:         %20%, DF: %21%\n\n"
                                 "   ***\n"
                                 "   *** The following statistics are computed on the log2 scale.\n"
                                 "   ***\n"
                                 "   ***   Eg: If the data points are (1,1), (2,2). The correlation will\n"
                                 "   ***       be computed on (log2(1), log2(1)), (log2(2), log2(2)))\n"
                                 "   ***\n\n"
                                 "   Correlation: %22%\n"
                                 "   Slope:       %23%\n"
                                 "   R2:          %24%\n"
                                 "   F-statistic: %25%\n"
                                 "   P-value:     %26%\n"
                                 "   SSM:         %27%, DF: %28%\n"
                                 "   SSE:         %29%, DF: %30%\n"
                                 "   SST:         %31%, DF: %32%\n";

            const auto n_lm = stats.linear(false);
            const auto l_lm = stats.linear(true);

            writer->open(file);
            writer->write((boost::format(summary) % data                         // 1
                                                  % stats.n_expT
                                                  % (samples.empty() ? "Experiment" : samples)
                                                  % stats.n_chrT
                                                  % units                        // 5
                                                  % stats.h.size()
                                                  % (ref.empty() ? units : ref)
                                                  % stats.ss.abund
                                                  % stats.ss.id
                                                  % detectHist(stats.h)           // 10
                                                  % n_lm.r
                                                  % n_lm.m
                                                  % n_lm.r2
                                                  % n_lm.f
                                                  % n_lm.p
                                                  % n_lm.ssm
                                                  % n_lm.ssm_df
                                                  % n_lm.sse
                                                  % n_lm.sse_df
                                                  % n_lm.sst
                                                  % n_lm.sst_df
                                                  % l_lm.r                        // 22
                                                  % l_lm.m
                                                  % l_lm.r2
                                                  % l_lm.f
                                                  % l_lm.p
                                                  % l_lm.ssm
                                                  % l_lm.ssm_df
                                                  % l_lm.sse
                                                  % l_lm.sse_df
                                                  % l_lm.sst
                                                  % l_lm.sst_df
                                                 ).str());
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

        /*
         * Provides a common framework to generate a CSV for all sequins
         */
        
        template <typename Writer> static void writeCSV(const std::vector<double> &x,
                                                        const std::vector<double> &y,
                                                        const std::vector<std::string> &z,
                                                        const FileName &file,
                                                        const std::string &xLabel,
                                                        const std::string &yLabel,
                                                        Writer writer)
        {
            writer->open(file);
            writer->write((boost::format("ID\t%1%\t%2%") % xLabel % yLabel).str());

            /*
             * Prefer to write results in sorted order
             */

            std::set<std::string> sorted(z.begin(), z.end());

            for (const auto &s : sorted)
            {
                const auto it = std::find(z.begin(), z.end(), s);
                const auto i  = std::distance(z.begin(), it);

                writer->write((boost::format("%1%\t%2%\t%3%") % z.at(i) % x.at(i) % y.at(i)).str());
            }

            writer->close();
        }

        /*
         * Provides a common framework for generating a R scatter plot
         */

        template <typename Stats, typename Writer> static void
                        scatter(const Stats &stats,
                                const std::string &title,
                                const std::string &prefix,
                                const AxisLabel &xLabel,
                                const AxisLabel &yLabel,
                                const AxisLabel &xLogLabel,
                                const AxisLabel &yLogLabel,
                                Writer writer,
                                bool shoudLog2 = true,
                                bool shouldCSV = true)
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
             * Generate an R script for data visualization
             */

            writer->open(prefix + "_plot.R");
            writer->write(RWriter::coverage(x, y, z, shoudLog2 ? xLogLabel : xLabel, shoudLog2 ? yLogLabel : yLabel, title, stats.s.abund));
            writer->close();

            /*
             * Generate CSV for each sequin
             */

            if (shouldCSV)
            {
                writeCSV(x, y, z, prefix + "_quin.csv", xLabel, yLabel, writer);                
            }
        }
    };
}

#endif
