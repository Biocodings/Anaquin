#ifndef STATS_HPP
#define STATS_HPP

#include <vector>
#include <iomanip>
#include <sstream>
#include "data/data.hpp"
#include <boost/format.hpp>
#include <ss/maths/maths.hpp>
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
        return (join + boost::algorithm::join(x, join + ", " + join) + join);
    }

    #define STRING(x) static_cast<std::string>(x)

    template <typename T> class SSamples
    {
        public:
        
            inline void add(const T &x)
            {
                _data.push_back(x);
            }
        
            inline std::size_t size() const { return _data.size(); }
        
            virtual operator std::string() const = 0;

        protected:
        
            std::vector<T> _data;
    };
    
    struct SReals : public SSamples<double>
    {
        inline double getMean() const { return SS::getMean(SSamples<double>::_data); }
        
        virtual operator std::string() const
        {
            const auto &data = SSamples<double>::_data;
            
            if (data.size() == 0)
            {
                return "-";
            }
            else if (data.size() > 1)
            {
                return (boost::format("%1$.2f \u00B1 %2$.2f") % SS::getMean(data) % SS::getSD(data)).str();
            }
            else
            {
                return (boost::format("%1$.2f") % data.front()).str();
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

    struct SPercent : public SReals
    {
        virtual operator std::string() const
        {
            return SReals::operator std::string() + "%";
        }
    };
    
    struct SProps : public SReals
    {
        // Emtpy Implementation
    };
    
    struct SCounts : public SSamples<Counts>
    {
        virtual operator std::string() const
        {
            const auto &data = SSamples<Counts>::_data;
            
            if (data.size() == 0)
            {
                return "-";
            }
            else if (data.size() > 1)
            {
                return (boost::format("%1% \u00B1 %2%") % SS::getMean(data) % SS::getSD(data)).str();
            }
            else
            {
                return (boost::format("%1%") % data.front()).str();
            }
        }
    };
}

#endif