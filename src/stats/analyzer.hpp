#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <map>
#include <memory>
#include <numeric>
#include "stats/limit.hpp"
#include <boost/format.hpp>
#include "stats/classify.hpp"
#include "writers/r_writer.hpp"
#include "writers/mock_writer.hpp"
#include <ss/regression/linear.hpp>

namespace Anaquin
{
    typedef std::map<BinID, Counts> BinCounts;

    template <typename T> std::string toNA(const T &t)
    {
        return isnan(t) ? "NA" : std::to_string(t);
    }
    
    inline std::string extractFile(const std::string &path)
    {
        auto r   = path;
        auto sep = r.find_last_of("\\/");
        
        if (sep != std::string::npos)
        {
            r = path.substr(sep + 1, path.size() - sep - 1);
        }

        return r;
    }
    
    template <typename T> Counts sum(const std::map<T, Counts> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), 0, [](Counts c, const std::pair<T, Counts>& p)
        {
            return c + p.second;
        });
    }

    struct Expression
    {
        ChromoID cID;
        
        GenericID id;
        
        Locus l;

        // The expression for signal
        FPKM fpkm;
    };

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

        // Total mapped to the synthetic chromosome
        Counts n_chrT = 0;

        // Total mapped to the endogenous
        Counts n_endo = 0;

        // Proportion of reads aligned to endogenous
        inline Percentage endoProp() const
        {
            return (n_chrT + n_endo) ? static_cast<double>(n_endo) / (n_chrT + n_endo) : NAN;
        }
        
        // Proportion of reads aligned to synthetic
        inline Percentage chrTProp() const
        {
            return dilution();
        }
        
        inline Percentage dilution() const
        {
            return (n_chrT + n_endo) ? static_cast<double>(n_chrT) / (n_chrT + n_endo) : NAN;
        }
    };

    struct AlignmentStats : public MappingStats
    {
        Counts unmapped = 0;

        template <typename T, typename F> void update(const T &t, F f)
        {
            if (!t.i)
            {
                if      (!t.mapped) { unmapped++; }
                else if (!f(t))     { n_chrT++;   }
                else                { n_endo++;   }
            }
        }

        template <typename T> void update(const T &t)
        {
            return update(t, [&](const T &t)
            {
                return t.id != ChrT;
            });
        }
    };

    /*
     * Represents something that is missing or undetected. It could be an exon, intron, isoform, gene etc.
     */
    
    struct Missing
    {
        Missing(const GenericID &id) : id(id) {}
        
        inline bool operator==(const Missing &m) const { return id == m.id; }
        inline bool operator< (const Missing &m) const { return id <  m.id; }

        const GenericID id;
    };
    
    struct CountPercent
    {
        CountPercent(Counts i, Counts n) : i(i), n(n) {}
        
        inline operator std::string() const
        {
            return (boost::format("%1% (%2%)") % i % n).str();
        }

        inline double percent() const
        {
            assert(i <= n);
            return static_cast<double>(i) / n;
        }
        
        // Relevant number of counts
        Counts i;
        
        // Total number of counts
        Counts n;        
    };

    struct UnknownAlignment
    {
        UnknownAlignment(const std::string &id, const Locus &l) : l(l) {}

        // Eg: HISEQ:132:C7F8BANXX:7:1116:11878:9591
        std::string id;
        
        // The position of the alignment
        Locus l;
    };
    
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
    
    /*
     * Represents an inflection limit, typically estimated by segmented piecewise linear regression.
     */
    
    struct InflectionLimit
    {
        // Name of the sequin
        std::string id;
        
        // Coefficient of determination before and after the break-point
        double lR2, rR2;
        
        // Slope before and after the break-point
        double lSl, rSl;
        
        // Intercept before and after the break-point
        double lInt, rInt;
        
        // The optimal breakpoint
        double b;
    };

    // Classify at the base-level by counting for non-overlapping regions
    template <typename I1, typename I2> void countBase(const I1 &r, const I2 &q, Confusion &m, SequinHist &c)
    {
        typedef typename I2::value_type Type;
        
        assert(!Locus::overlap(r));
        
        const auto merged = Locus::merge<Type, Locus>(q);
        
        for (const auto &l : merged)
        {
            m.nq() += l.length();
            m.tp() += countOverlaps(r, l, c);
            m.fp()  = m.nq() - m.tp();
            
            // Make sure we don't run into negative
            assert(m.nq() >= m.tp());
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

        // Return the x-values and y-values after filtering
        void data(std::vector<double> &x, std::vector<double> &y, bool shouldLog, std::vector<std::string> *ids = nullptr) const;

        // Compute the inflection limit. By default, this function assumes log-transformation.
        InflectionLimit inflect(bool shouldLog = true) const;
        
        // Compute a simple linear regression model. By default, this function assumes log-transformation.
        LinearModel linear(bool shouldLog = true) const;
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

    // Forward delcaration. Any analyzer doesn't need it need not link against it.
    class Experiment;
    
    struct AnalyzerOptions : public WriterOptions
    {
        // Experimental meta-data
        std::shared_ptr<Experiment> exp;

        // Only process this number of filters
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

    template <typename T> class Accumulator
    {
        public:
        
            typedef std::string Key;
        
            struct Deviation
            {
                // First moment
                T mean;
            
                // Standard deviation
                T sd;

                inline std::string operator()() const
                {
                    return (boost::format("%1% \u00B1 %2%") % mean % sd).str();
                }
            };
        
            void add(const Key &key, double value)
            {
                _data[key].push_back(value);
            }
        
            void add(const Key &key, const Limit &l)
            {
                _limits[key].push_back(l);
            }
        
            Deviation value(const Key &key) const
            {
                Deviation d;
            
                d.sd   = SS::sd(_data.at(key));
                d.mean = SS::mean(_data.at(key));
            
                return d;
            }
        
            const Limit & limits(const Key &key) const
            {
                const Limit *min = nullptr;
            
                for (const auto &limit : _limits.at(key))
                {
                    if (!min || limit.abund < min->abund)
                    {
                        min = &limit;
                    }
                }
            
                return *min;
            }
        
        private:
        
            // Used for comparing limit of detection
            std::map<Key, std::vector<Limit>> _limits;
        
            // Used for regular mappings
            std::map<Key, std::vector<double>> _data;
    };

    struct AnalyzeReporter
    {
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

    };
}

#endif
