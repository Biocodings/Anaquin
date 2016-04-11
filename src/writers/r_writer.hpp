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
    // Number of elements in the histogram with at least an entry
    template <typename T> Counts detect(const std::map<T, Counts> &m)
    {
        return std::count_if(m.begin(), m.end(), [&](const std::pair<T, long> &i)
        {
            return i.second ? 1 : 0;
        });
    }
    
    class MappingStats;
    
    struct StatsWriter
    {
        static Scripts writeCSV(const LinearStats &stats,
                                const Label &xLabel = "EAbund",
                                const Label &yLabel = "MAbund",
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
            ss << ((boost::format("Sequin\t%1%\t%2%\n") % xLabel % yLabel).str());

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
         * -------------------- Linear Statistics (with inflection) --------------------
         */

        typedef std::map<std::string, Counts> Hist;
        
        static Scripts inflectSummary();
        
        static Scripts inflectSummary(const FileName &,
                                      const FileName &,
                                      const SInflectStats &,
                                      const Units &);

        static Scripts inflectSummary(const FileName &,
                                      const FileName &,
                                      const std::vector<FileName>     &,
                                      const std::vector<Hist>         &,
                                      const std::vector<MappingStats> &,
                                      const std::vector<LinearStats>  &,
                                      const Units &);

        static Scripts inflectSummary(const FileName &,
                                      const FileName &,
                                      const FileName &,
                                      const Hist     &,
                                      const MappingStats &,
                                      const LinearStats  &,
                                      const Units &units);

        static Scripts inflectSummary(const FileName &,
                                      const FileName &,
                                      const Hist     &,
                                      const LinearStats  &,
                                      const Units &units);
        
        /*
         * -------------------- Linear Statistics --------------------
         */
        
        static Scripts linearSummary(const FileName &,
                                     const FileName &,
                                     const LinearStats &,
                                     const Hist &hist);

        template <typename Stats_1, typename Stats_2, typename Stats> static Scripts linear(const FileName &f,
                                                                                            const FileName &d1,
                                                                                            const FileName &d2,
                                                                                            const Stats_1  &s1,
                                                                                            const Stats_2  &s2,
                                                                                            const Stats    &s,
                                                                                            const ChrID    &cID,
                                                                                            const Units    &units,
                                                                                            const Units    &ref = "",
                                                                                            const Label    &samples = "")
        {
            const auto summary = "Summary for input: %1% and %2%\n\n"
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
            
            const auto n_lm = s.data.at(cID).linear(false);
            const auto l_lm = s.data.at(cID).linear(true);
            
            return (boost::format(summary) % d1                          // 1
                                           % d2
                                           % (samples.empty() ? "Genome" : samples)
                                           % s1.n_endo
                                           % units
                                           % s1.n_chrT
                                           % s2.n_endo
                                           % s2.n_chrT
                                           % s1.data.at(ChrT).h.size()
                                           % (ref.empty() ? units : ref) // 10
                                           % s2.data.at(ChrT).h.size()
                                           % detect(s1.data.at(ChrT).h)  // 12
                                           % detect(s2.data.at(ChrT).h)  // 13
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
                    ).str();
        }
    };
    
    struct RWriter
    {
        static Scripts createScript(const FileName &name, const Scripts &scripts);
        
        static Scripts createMA(const FileName &, const std::string &);
    };
}

#endif