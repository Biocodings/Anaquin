#ifndef GI_R_WRITER_HPP
#define GI_R_WRITER_HPP

#include <string>
#include <sstream>
#include <numeric>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

extern std::string linearR();

namespace Spike
{
    struct RWriter
    {
        template <typename I1, typename I2> static std::string write(const I1 &x, const I1 &y, const I2 &z)
        {
            using boost::algorithm::join;
            using boost::adaptors::transformed;

            std::stringstream ss;
            ss << linearR();

            std::cout << "Generating a linear model" << std::endl;
            
            return (boost::format(ss.str()) %
                        join(x | transformed(static_cast<std::string(*)(double)>(std::to_string)), ", ") %
                        join(y | transformed(static_cast<std::string(*)(double)>(std::to_string)), ", ") %
                       (boost::format("'%1%'") % boost::algorithm::join(z, "','")).str()).str();
        }
    };
}

#endif