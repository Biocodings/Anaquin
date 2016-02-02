#ifndef R_WRITER_HPP
#define R_WRITER_HPP

#include <math.h>
#include <sstream>
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
        /*
         * -------------------- Linear Statistics (with inflection) --------------------
         */

        static Scripts inflectSummary();
        static Scripts inflectSummary(const FileName &ref, const SInflectStats &stats);
        static Scripts inflectSummary(const FileName &ref,
                                      const std::vector<FileName> &,
                                      const std::vector<MappingStats> &,
                                      const std::vector<LinearStats> &,
                                      const Units &units);

        /*
         * -------------------- Linear Statistics --------------------
         */
        
        template <typename Stats> static Scripts linear(const FileName &src,
                                                        const Stats &stats,
                                                        const ChromoID &cID,
                                                        const Units &units,
                                                        const Units &ref = "")
        {
            const auto summary = "Summary for file: %1%\n\n"
                                 "   Experiment:  %2% %4%\n"
                                 "   Synthetic:   %3% %4%\n\n"
                                 "   Reference:   %5% %6%\n"
                                 "   Detected:    %9% %6%\n\n"
                                 "   ***\n"
                                 "   *** Detection Limits\n"
                                 "   ***\n\n"
                                 "   Absolute:    %7% (attomol/ul) (%8%)\n\n"
                                 "   ***\n"
                                 "   *** Statistics for linear regression\n"
                                 "   ***\n\n"
                                 "   Correlation: %10%\n"
                                 "   Slope:       %11%\n"
                                 "   R2:          %12%\n"
                                 "   F-statistic: %13%\n"
                                 "   P-value:     %14%\n"
                                 "   SSM:         %15%, DF: %16%\n"
                                 "   SSE:         %17%, DF: %18%\n"
                                 "   SST:         %19%, DF: %20%\n\n"
                                 "   ***\n"
                                 "   *** Statistics for linear regression (log2 scale)\n"
                                 "   ***\n\n"
                                 "   Correlation: %21%\n"
                                 "   Slope:       %22%\n"
                                 "   R2:          %23%\n"
                                 "   F-statistic: %24%\n"
                                 "   P-value:     %25%\n"
                                 "   SSM:         %26%, DF: %27%\n"
                                 "   SSE:         %28%, DF: %29%\n"
                                 "   SST:         %30%, DF: %31%\n";
            
            const auto n_lm = stats.data.at(cID).linear(false);
            const auto l_lm = stats.data.at(cID).linear(true);
            
            return (boost::format(summary) % src                          // 1
                                           % stats.n_endo
                                           % stats.n_chrT
                                           % units
                                           % stats.hist.size()            // 5
                                           % (ref.empty() ? units : ref)
                                           % stats.limit.abund
                                           % stats.limit.id
                                           % detect(stats.hist)
                                           % n_lm.r                       // 10
                                           % n_lm.m
                                           % n_lm.R2
                                           % n_lm.F
                                           % n_lm.p
                                           % n_lm.SSM
                                           % n_lm.SSM_D
                                           % n_lm.SSE
                                           % n_lm.SSE_D
                                           % n_lm.SST
                                           % n_lm.SST_D
                                           % l_lm.r
                                           % l_lm.m                       // 22
                                           % l_lm.R2
                                           % l_lm.F
                                           % l_lm.p
                                           % l_lm.SSM
                                           % l_lm.SSM_D
                                           % l_lm.SSE
                                           % l_lm.SSE_D
                                           % l_lm.SST
                                           % l_lm.SST_D                  // 31
                    ).str();
        }
        
        /*
         * -------------------- Linear Statistics --------------------
         */
        
        template <typename Stats_1, typename Stats_2, typename Stats> static Scripts linear(const FileName &f,
                                                                                            const FileName &d1,
                                                                                            const FileName &d2,
                                                                                            const Stats_1  &s1,
                                                                                            const Stats_2  &s2,
                                                                                            const Stats    &s,
                                                                                            const ChromoID &cID,
                                                                                            const Units    &units,
                                                                                            const Units &ref = "",
                                                                                            const Label &samples = "")
        {
            const auto summary = "Summary for file: %1% and %2%\n\n"
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
        /*
         * -------------------- ROC Plot --------------------
         */

        // Create an ROC script for transcriptome analysis
        static Scripts createROC_T(const std::vector<FeatureID> &, const std::vector<double> &, const std::string &);

        // Create an ROC script for variant analysis
        static Scripts createROC_V(const Path &, const FileName &, const FileName &);

        /*
         * -------------------- MA Plot --------------------
         */
        
        // Create an MA plot and link to the count table
        static Scripts createMA(const Path &, const FileName &, const std::string &);

        /*
         * -------------------- LODR Plot --------------------
         */
        
        // Create a LODR script for transcriptome analysis
        static Scripts createLODR_T(const Path &, const FileName &);

        // Create a LODR script for variant analysis
        static Scripts createLODR_V(const Path &, const FileName &);

        /*
         * --------------------- Splice Plot ---------------------
         */

        static Scripts createSplice(const Path &, const FileName &);

        /*
         * -------------------- Scatter Plot --------------------
         */

        static Scripts scatterPool(const Path &, const FileName &);

        template <typename Stats> static Scripts scatter(const Stats &stats,
                                                         const ChromoID &cID,
                                                         const std::string &title,
                                                         const std::string &prefix,
                                                         const AxisLabel &xLabel,
                                                         const AxisLabel &yLabel,
                                                         const AxisLabel &xLogLabel,
                                                         const AxisLabel &yLogLabel,
                                                         bool  shoudLog2 = true)
        {
            std::vector<double> x, y;
            std::vector<std::string> z;
            
            /*
             * Ignore any invalid value...
             */
            
            for (const auto &p : stats.data.at(cID))
            {
                if (!isnan(p.second.x) && !isnan(p.second.y))
                {
                    z.push_back(p.first);
                    x.push_back(p.second.x);
                    y.push_back(p.second.y);
                }
            }
            
            //return RWriter::scatter(z, x, y, shoudLog2 ? xLogLabel : xLabel, shoudLog2 ? yLogLabel : yLabel, title, stats.limit.abund);
            return RWriter::scatter(z, x, y, shoudLog2 ? xLogLabel : xLabel, shoudLog2 ? yLogLabel : yLabel, title);
        }

        template <typename T> static Scripts scatter(const T &t,
                                                     const std::string &title,
                                                     const AxisLabel &xLabel,
                                                     const AxisLabel &yLabel,
                                                     const AxisLabel &xLogLabel,
                                                     const AxisLabel &yLogLabel,
                                                     bool  shoudLog2 = true)
        {
            std::vector<double> x, y;
            std::vector<std::string> z;

            for (const auto &p : t)
            {
                if (!isnan(p.second.x) && !isnan(p.second.y))
                {
                    z.push_back(p.first);
                    x.push_back(p.second.x);
                    y.push_back(p.second.y);
                }
            }
            
            return RWriter::scatter(z, x, y, shoudLog2 ? xLogLabel : xLabel, shoudLog2 ? yLogLabel : yLabel, title);
        }
        
        static Scripts scatter(const std::vector<SequinID> &seqs,
                               const std::vector<double>   &x,
                               const std::vector<double>   &y,
                               const AxisLabel &xLabel,
                               const AxisLabel &yLabel,
                               const AxisLabel &title);
    };
}

#endif