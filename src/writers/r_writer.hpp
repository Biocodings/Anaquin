#ifndef R_WRITER_HPP
#define R_WRITER_HPP

#include <map>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
#include <numeric>
#include <iomanip>
#include "data/types.hpp"
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

// Defined in main.cpp
extern std::string date();

// Defined in main.cpp
extern std::string __full_command__;

// Defined in resources.cpp
extern std::string PlotScatter();

// Defined in resources.cpp
extern std::string PlotMA();

// Defined in resources.cpp
extern std::string PlotROC();

// Defined in resources.cpp
extern std::string PlotLODR();

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
    
    struct StatsWriter
    {
        /*
         * -------------------- Linear Statistics (with inflection) --------------------
         */

        template <typename Stats> static std::string linearInflect(const FileName &src,
                                                                   const Stats &stats,
                                                                   const ChromoID &cID,
                                                                   const Units &units,
                                                                   const Units &ref = "")
        {
            const auto summary = "Summary for dataset: %1%\n\n"
                                 "   Experiment:  %2% %4%\n"
                                 "   Synthetic:   %3% %4%\n\n"
                                 "   Reference:   %5% %6%\n"
                                 "   Detected:    %9% %6%\n\n"
                                 "   ***\n"
                                 "   *** Detection Limits\n"
                                 "   ***\n\n"
                                 ""
                                 "   Break: %7% (%8%)\n\n"
                                 "   Left:  %9% + %10%x (R2 = %11%)\n"
                                 "   Right: %12% + %13%x (R2 = %14%)\n\n"
                                 "   ***\n"
                                 "   *** Statistics for linear regression\n"
                                 "   ***\n\n"
                                 "   Correlation: %15%\n"
                                 "   Slope:       %16%\n"
                                 "   R2:          %17%\n"
                                 "   F-statistic: %18%\n"
                                 "   P-value:     %19%\n"
                                 "   SSM:         %20%, DF: %21%\n"
                                 "   SSE:         %22%, DF: %23%\n"
                                 "   SST:         %24%, DF: %25%\n\n"
                                 "   ***\n"
                                 "   *** Statistics for linear regression (log2 scale)\n"
                                 "   ***\n\n"
                                 "   Correlation: %26%\n"
                                 "   Slope:       %27%\n"
                                 "   R2:          %28%\n"
                                 "   F-statistic: %29%\n"
                                 "   P-value:     %30%\n"
                                 "   SSM:         %31%, DF: %32%\n"
                                 "   SSE:         %33%, DF: %34%\n"
                                 "   SST:         %35%, DF: %36%\n";
            
            const auto n_lm = stats.data.at(cID).linear(false);
            const auto l_lm = stats.data.at(cID).linear(true);
            
            // Calcluate the inflect point after log-transformation
            const auto inflect = stats.data.at(cID).inflect(true);
            
            // Remember the break-point is on the log-scale, we'll need to convert it back
            const auto b = pow(2, inflect.b);
            
            return (boost::format(summary) % src                          // 1
                                           % stats.n_endo
                                           % stats.n_chrT
                                           % units
                                           % stats.hist.size()            // 5
                                           % (ref.empty() ? units : ref)  // 6
                                           % b
                                           % inflect.id
                                           % inflect.lInt                 // 9
                                           % inflect.lSl                  // 10
                                           % inflect.lR2                  // 11
                                           % inflect.rInt                 // 12
                                           % inflect.rSl                  // 13
                                           % inflect.rR2                  // 14
                                           % n_lm.r                       // 15
                                           % n_lm.m                       // 16
                                           % n_lm.r2                      // 17
                                           % n_lm.f                       // 18
                                           % n_lm.p                       // 19
                                           % n_lm.ssm                     // 20
                                           % n_lm.ssm_df                  // 21
                                           % n_lm.sse                     // 22
                                           % n_lm.sse_df                  // 23
                                           % n_lm.sst                     // 24
                                           % n_lm.sst_df                  // 25
                                           % l_lm.r                       // 26
                                           % l_lm.m                       // 27
                                           % l_lm.r2                      // 28
                                           % l_lm.f                       // 29
                                           % l_lm.p                       // 30
                                           % l_lm.ssm                     // 31
                                           % l_lm.ssm_df                  // 32
                                           % l_lm.sse                     // 33
                                           % l_lm.sse_df                  // 34
                                           % l_lm.sst                     // 35
                                           % l_lm.sst_df                  // 36
                    ).str();
        }
        
        
        
        /*
         * -------------------- Linear Statistics --------------------
         */
        
        template <typename Stats> static std::string linear(const FileName &src,
                                                            const Stats &stats,
                                                            const ChromoID &cID,
                                                            const Units &units,
                                                            const Units &ref = "")
        {
            const auto summary = "Summary for dataset: %1%\n\n"
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
                                           % n_lm.r2
                                           % n_lm.f
                                           % n_lm.p
                                           % n_lm.ssm
                                           % n_lm.ssm_df
                                           % n_lm.sse
                                           % n_lm.sse_df
                                           % n_lm.sst
                                           % n_lm.sst_df
                                           % l_lm.r
                                           % l_lm.m                       // 22
                                           % l_lm.r2
                                           % l_lm.f
                                           % l_lm.p
                                           % l_lm.ssm
                                           % l_lm.ssm_df
                                           % l_lm.sse
                                           % l_lm.sse_df
                                           % l_lm.sst
                                           % l_lm.sst_df                  // 31
                    ).str();
        }
        
        /*
         * -------------------- Linear Statistics --------------------
         */
        
        template <typename Stats_1, typename Stats_2, typename Stats> static std::string linear(const FileName &f,
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
    
    class CountTable;
    
    struct RWriter
    {
        static std::string d2str(double x)
        {
            std::ostringstream out;
            out << std::setprecision(6) << x;
            return out.str();
        }
        
        static std::string x2str(unsigned x)
        {
            return std::to_string(x);
        }
        
        template <typename T> static std::string concat(const std::vector<T> &x, std::string (*f)(T) = d2str)
        {
            return boost::algorithm::join(x | boost::adaptors::transformed(static_cast<std::string(*)(T)>(f)), ", ");
        }

        static std::string concat(const std::vector<std::string> &x)
        {
            return ("\'" + boost::algorithm::join(x, "\',\'") + "\'");
        }

        /*
         * -------------------- ROC Plot --------------------
         */
        
        static Scripts createROC(const std::vector<std::string> &, const std::vector<double> &, const std::string &);

        /*
         * -------------------- MA Plot --------------------
         */
        
        // Create an MA plot and link to the count table
        static Scripts createMA(const std::string &working, const FileName &, const std::string &);

        /*
         * -------------------- LODR Plot --------------------
         */
        
        static Scripts createLODR(const std::string &, const FileName &, const FileName &, const std::string &);

        /*
         * -------------------- Scatter Plot --------------------
         */
        
        template <typename Stats> static Scripts scatter(const Stats &stats,
                                                         const ChromoID &cID,
                                                         const std::string &title,
                                                         const std::string &prefix,
                                                         const AxisLabel &xLabel,
                                                         const AxisLabel &yLabel,
                                                         const AxisLabel &xLogLabel,
                                                         const AxisLabel &yLogLabel,
                                                         bool shoudLog2 = true)
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
            
            return RWriter::scatter(z, x, y, shoudLog2 ? xLogLabel : xLabel, shoudLog2 ? yLogLabel : yLabel, title, stats.limit.abund);
        }
        
        template <typename T> static Scripts scatter(const std::vector<SequinID> &seqs,
                                                     const std::vector<T>        &x,
                                                     const std::vector<T>        &y,
                                                     const std::string           &xLabel,
                                                     const std::string           &yLabel,
                                                     const std::string           &title,
                                                     T s)
        {
            assert(!xLabel.empty() && !yLabel.empty());

            std::stringstream ss;
            ss << PlotScatter();

            return (boost::format(ss.str()) % date()
                                            % __full_command__
                                            % ("\'" + boost::algorithm::join(seqs, "\',\'") + "\'")
                                            % RWriter::concat(x)
                                            % RWriter::concat(y)
                                            % xLabel
                                            % yLabel).str();
        }
    };
}

#endif