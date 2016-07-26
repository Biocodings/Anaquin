#ifndef R_WRITER_HPP
#define R_WRITER_HPP

#include <set>
#include <math.h>
#include <numeric>
#include "stats/linear.hpp"
#include <boost/format.hpp>

// Defined in main.cpp
extern std::string date();

// Defined in main.cpp
extern std::string __full_command__;

namespace Anaquin
{
    class MappingStats;
    
    struct StatsWriter
    {
        typedef std::map<std::string, Counts> Hist;

        static SInflectStats multiInfect(const std::vector<FileName>     &,
                                         const std::vector<MappingStats> &,
                                         const std::vector<LinearStats>  &);
        
        static Scripts writeCSV(const LinearStats &stats,
                                const Label &xLabel = "Expected",
                                const Label &yLabel = "Measured",
                                bool shouldLog = false)
        {
            const auto d = stats.data(false);
            return StatsWriter::writeCSV(d.x, d.y, d.ids, xLabel, yLabel, shouldLog);
        }

        static Scripts writeCSV(const std::vector<double> &x,
                                const std::vector<double> &y,
                                const std::vector<SequinID> &ids,
                                const Label &xLabel,
                                const Label &yLabel,
                                bool shouldLog = false)
        {
            std::stringstream ss;
            ss << ((boost::format("ID\t%1%\t%2%\n") % xLabel % yLabel).str());

            std::set<SequinID> sorted(ids.begin(), ids.end());
            
            for (const auto &s : sorted)
            {
                const auto it = std::find(ids.begin(), ids.end(), s);
                const auto i  = std::distance(ids.begin(), it);
                
                ss << ((boost::format("%1%\t%2%\t%3%\n") % ids.at(i)
                                                         % (shouldLog ? log2(x.at(i)) : x.at(i))
                                                         % (shouldLog ? log2(y.at(i)) : y.at(i))).str());
            }

            return ss.str();
        }

        /*
         * -------------------- Linear Statistics --------------------
         */
        
        static Scripts linearSummary(const FileName &,
                                     const FileName &,
                                     const LinearStats &,
                                     const MappingStats &,
                                     const Hist  &,
                                     const Units &units);
    };
    
    struct RWriter
    {
        static Scripts createVROC(const FileName &, const std::string &);
        
        static Scripts createSensitivity(const FileName    &,
                                         const std::string &,
                                         const std::string &,
                                         const std::string &,
                                         const std::string &,
                                         const std::string &,
                                         bool showLOQ);
        
        static Scripts createMultiScatter(const FileName  &,
                                          const std::string &,
                                          const std::string &,
                                          const std::string &,
                                          const std::string &,
                                          const std::string &,
                                          const std::string &,
                                          bool showLOQ,
                                          bool shouldLog);

        static Scripts createFold(const FileName    &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  bool shouldLog);
        
        static Scripts createScatterNeedLog(const FileName  &,
                                            const std::string &,
                                            const std::string &,
                                            const std::string &,
                                            const std::string &,
                                            const std::string &,
                                            const std::string &,                                            
                                            bool showLOQ);

        static Scripts createScatterNoLog(const FileName    &,
                                          const std::string &,
                                          const std::string &,
                                          const std::string &,
                                          const std::string &,
                                          const std::string &,
                                          const std::string &,
                                          bool showLOQ);
        
        static Scripts createScript(const FileName &, const Scripts &);
        static Scripts createScript(const FileName &, const Scripts &, const std::string &);
    };
}

#endif