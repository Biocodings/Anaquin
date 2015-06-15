#ifndef GI_R_WRITER_HPP
#define GI_R_WRITER_HPP

#include <string>
#include <vector>
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
        template <typename T> static std::string write
                    (const std::vector<T> &x,
                     const std::vector<T> &y,
                     const std::vector<SequinID> &z)
        {
            std::vector<ColorID> c(z.size(), "black");

            // Apply the default drawing color
            return RWriter::write(x, y, z, c);
        }
        
        template <typename T> static std::string write
                    (const std::vector<T> &x,
                     const std::vector<T> &y,
                     const std::vector<SequinID> &z,
                     const std::vector<ColorID>  &c)
        {
            using boost::algorithm::join;
            using boost::adaptors::transformed;

            std::stringstream ss;
            ss << linearR();

            const auto xs = join(x | transformed(static_cast<std::string(*)(double)>(std::to_string)), ", ");
            const auto ys = join(y | transformed(static_cast<std::string(*)(double)>(std::to_string)), ", ");
            const auto zs = (boost::format("'%1%'") % boost::algorithm::join(z, "','")).str();
            const auto cs = (boost::format("'%1%'") % boost::algorithm::join(c, "','")).str();

            return (boost::format(ss.str()) % xs % ys % zs % cs).str();
        }
    };
}

#endif