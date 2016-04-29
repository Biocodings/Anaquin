#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include "stats/classify.hpp"
#include "writers/r_writer.hpp"
#include "writers/mock_writer.hpp"

namespace Anaquin
{
    template <typename T> std::string toString(const T &x, unsigned n = 2)
    {
        std::ostringstream out;
        out << std::fixed << std::setprecision(n) << x;
        return out.str();
    }

    template <typename T> std::string n2str(const T &x)
    {
        return isnan(x) ? "NA" : toString(x);
    }
    
    template <typename T> std::string p2str(const T &x)
    {
        return isnan(x) ? "NA" : toString(x, 310);
    }
    
    template <typename T> Counts count(const std::map<T, Counts> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), 0, [](Counts c, const std::pair<T, Counts>& p)
        {
            return c + (p.second ? 1 : 0);
        });
    }
    
    template <typename T> Counts sum(const std::map<T, Counts> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), 0, [](Counts c, const std::pair<T, Counts>& p)
        {
            return c + p.second;
        });
    }

    struct Analyzer
    {
        // Empty Implementation
    };
    
    typedef std::map<SequinID, Counts> Hist;
    
    struct SequinStats
    {
        // Distribution of counts across sequins
        Hist hist;
    };

    struct MappingStats
    {
        inline Proportion endoProp() const
        {
            return (n_chrT + n_geno) ? static_cast<Proportion>(n_geno) / (n_chrT + n_geno) : NAN;
        }

        inline Proportion chrTProp() const
        {
            return dilution();
        }

        inline Proportion dilution() const
        {
            return (n_chrT + n_geno) ? static_cast<Proportion>(n_chrT) / (n_chrT + n_geno) : NAN;
        }

        // Total mapped to the synthetic chromosome
        Counts n_chrT = 0;

        // Total mapped to the endogenous
        Counts n_geno = 0;
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
                else                { n_geno++;   }
            }
        }

        template <typename T> void update(const T &t)
        {
            return update(t, [&](const T &t)
            {
                return t.cID != ChrT;
            });
        }
    };

    struct FusionStats : public MappingStats
    {
        // Number of fusions spanning across the genome and the synthetic chromosome
        Counts chrT_endo = 0;
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
        Path work;
        
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
        // Reference annotation for the synthetic
        FileName rChrT;
        
        // Reference annotation for the endogenous
        FileName rGeno;
    };
    
    struct IndexOptions : public AnalyzerOptions
    {
        FileName index;
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
    
    struct ReportOptions : public AnalyzerOptions
    {
        FileName mix;
        FileName index;
    };

    struct DoubleMixtureOptions : public AnalyzerOptions
    {
        Mixture mix_1 = Mix_1;
        Mixture mix_2 = Mix_2;
    };
    
    struct CalledFusion
    {
        inline operator Locus() const
        {
            return Locus(l1, l2);
        }
        
        // Chromosome for the left and right
        ChrID cID_1, cID_2;
        
        // Strand for the left and right
        Strand s1, s2;
        
        // Where the fusion occurs
        Base l1, l2;
        
        // Number of reads that span the fusion (junction reads)
        Reads reads;
    };
}

#endif
