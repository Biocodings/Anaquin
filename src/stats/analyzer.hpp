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
            return nSeqs + nEndo + nNA;
        }

        inline Proportion pNA() const
        {
            return total() ? static_cast<Proportion>(100.0 * nNA) / total() : NAN;
        }
        
        inline Proportion pEndo() const
        {
            return total() ? static_cast<Proportion>(100.0 * nEndo) / total() : NAN;
        }

        inline Proportion pSyn() const
        {
            return total() ? static_cast<Proportion>(100.0 * nSeqs) / total() : NAN;
        }

        inline Proportion dilution() const
        {
            return (nSeqs + nEndo) ? static_cast<Proportion>(nSeqs) / (nSeqs + nEndo) : NAN;
        }

        Counts nNA   = 0;
        Counts nEndo = 0;
        Counts nSeqs = 0;
    };

    struct AlignmentStats : public MappingStats
    {
        template <typename T, typename F> void update(const T &t, F f)
        {
            if      (!t.mapped) { nNA++;   }
            else if (f(t.cID))  { nSeqs++; }
            else                { nEndo++; }
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

        std::shared_ptr<Writer<>> writer = std::shared_ptr<Writer<>>(new MockWriter());
        std::shared_ptr<Writer<>> logger = std::shared_ptr<Writer<>>(new MockWriter());
        std::shared_ptr<Writer<>> output = std::shared_ptr<Writer<>>(new MockWriter());

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

    struct FuzzyOptions : public AnalyzerOptions
    {
        double fuzzy;
    };
}

#endif
