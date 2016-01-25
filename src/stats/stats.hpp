#ifndef STATS_HPP
#define STATS_HPP

#include <vector>
#include <iomanip>
#include <sstream>
#include <ss/stats.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

namespace Anaquin
{
    static inline std::string d2str(double x)
    {
        std::ostringstream out;
        out << std::setprecision(6) << x;
        return out.str();
    }
    
    static inline std::string x2str(unsigned x)
    {
        return std::to_string(x);
    }

    template <typename T> static std::string concat(const std::vector<T> &x, std::string (*f)(T) = d2str)
    {
        return boost::algorithm::join(x | boost::adaptors::transformed(static_cast<std::string(*)(T)>(f)), ", ");
    }
    
    static std::string concat(const std::vector<std::string> &x, const std::string &join = "\'")
    {
        return (join + boost::algorithm::join(x, join + "," + join) + join);
    }

    #define STRING(x) static_cast<std::string>(x)

    template <typename T> class SSamples
    {
        public:
        
            void add(const T &x)
            {
                _data.push_back(x);
            }
        
            virtual operator std::string() const = 0;
        
        protected:
        
            std::vector<T> _data;
    };
    
    struct SReals : public SSamples<double>
    {
        virtual operator std::string() const
        {
            if (SSamples<double>::_data.size() == 0)
            {
                return "-";
            }
            else if (SSamples<double>::_data.size() > 1)
            {
                return (boost::format("%1% \u00B1 %2%") % SS::mean(SSamples<double>::_data)
                                                        % SS::sd(SSamples<double>::_data)).str();
            }
            else
            {
                return (boost::format("%1%") % SSamples<double>::_data.front()).str();
            }
        }
    };
    
    struct SStrings : public SSamples<std::string>
    {
        virtual operator std::string() const
        {
            return concat(SSamples<std::string>::_data, "");
        }
    };

    struct SProps : public SReals
    {
        virtual operator std::string() const
        {
            return SReals::operator std::string() + "%";
        }
    };
    
    typedef SReals SCounts;
}

#endif