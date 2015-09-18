#ifndef GI_R_WRITER_HPP
#define GI_R_WRITER_HPP

#include <string>
#include <vector>
#include <sstream>
#include <numeric>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

// Defined in main.cpp
extern std::string date();

// Defined in main.cpp
extern std::string __full_command__;

extern std::string RScriptCoverage();

namespace Anaquin
{
    struct RWriter
    {
        // Generate a R script for plotting expected and measured coverage
        template <typename T> static std::string coverage
                        (const std::vector<T> &x,
                         const std::vector<T> &y,
                         const std::vector<SequinID> &z,
                         const std::string &xLabel,
                         const std::string &yLabel,
                         T s)
        {
            assert(!xLabel.empty() && !yLabel.empty());

            using boost::algorithm::join;
            using boost::adaptors::transformed;

            std::stringstream ss;
            ss << RScriptCoverage();

            const auto xs = join(x | transformed(static_cast<std::string(*)(double)>(std::to_string)), ", ");
            const auto ys = join(y | transformed(static_cast<std::string(*)(double)>(std::to_string)), ", ");
            const auto zs = (boost::format("'%1%'") % boost::algorithm::join(z, "','")).str();

            return (boost::format(ss.str()) % date()
                                            % __full_command__
                                            % xs
                                            % ys
                                            % zs
                                            % xLabel
                                            % yLabel).str();
        }
    };
}

#endif