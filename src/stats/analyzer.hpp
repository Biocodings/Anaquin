#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include "data/standard.hpp"
#include "stats/classify.hpp"
#include "writers/r_writer.hpp"
#include "writers/mock_writer.hpp"

// Defined in main.cpp
extern bool __showInfo__;

namespace Anaquin
{
    struct Analyzer
    {
        // Empty Implementation
    };
    
    struct Sample
    {
        Path path;
        FileName p1, p2;
    };
    
    typedef std::map<SequinID, Counts> Hist;
    
    struct HistStats
    {
        // Distribution of counts within sampling regions
        Hist hist;
    };

    struct AnalyzerStats
    {
        // Empty Implementation
    };
    
    struct LimitStats
    {
        // Absolute detection limit
        Limit limit;
    };

    struct MappingStats
    {
        inline Counts total() const
        {
            return nSyn + nGen + nNA;
        }

        inline Proportion propNA() const
        {
            return total() ? static_cast<Proportion>(100.0 * nNA) / total() : NAN;
        }
        
        inline Proportion propGen() const
        {
            return total() ? static_cast<Proportion>(100.0 * nGen) / total() : NAN;
        }

        inline Proportion propSyn() const
        {
            return total() ? static_cast<Proportion>(100.0 * nSyn) / total() : NAN;
        }

        inline Proportion dilution() const
        {
            return (nSyn + nGen) ? static_cast<Proportion>(nSyn) / (nSyn + nGen) : NAN;
        }

        Counts nNA  = 0;
        Counts nGen = 0;
        Counts nSyn = 0;
    };

    struct AlignmentStats : public MappingStats
    {
        template <typename T, typename F> void update(const T &t, F f)
        {
            if      (!t.mapped) { nNA++;  }
            else if (f(t.cID))  { nSyn++; }
            else                { nGen++; }
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
        
        inline void generate(const FileName &f) const
        {
            info("Generating: " + f);
        }
        
        inline void info(const std::string &s) const
        {
            logInfo(s);
            
            if (__showInfo__)
            {
                output->write("[INFO]: " + s);
            }
        }
        
        inline void logInfo(const std::string &s) const
        {
            logger->write("[INFO]: " + s);
        }
        
        inline void logWarn(const std::string &s) const
        {
            logger->write("[WARN]: " + s);
        }

        inline void logWait(const std::string &s) const
        {
            logger->write("[WAIT]: " + s);
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
        // Empty Implementation
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
}

#endif
