#ifndef GI_R_WRITER_HPP
#define GI_R_WRITER_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

namespace Spike
{
    struct RWriter
    {
        template <typename Iter> static std::string write(Iter &x, Iter &y)
        {
            using boost::algorithm::join;
            using boost::adaptors::transformed;
            
            std::stringstream ss;
            ss << std::ifstream("scripts/exp_t.R").rdbuf();

            const auto xs = join(x | transformed(static_cast<std::string(*)(double)>(std::to_string)), ", " );
            const auto ys = join(y | transformed(static_cast<std::string(*)(double)>(std::to_string)), ", " );
    
            return (boost::format(ss.str()) % xs % ys).str();
        }
    };
}

#endif