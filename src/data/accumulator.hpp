#ifndef ACCUMULATOR_HPP
#define ACCUMULATOR_HPP

#include <map>
#include <ss/stats.hpp>
#include "stats/limit.hpp"
#include <boost/format.hpp>

namespace Anaquin
{
    class Accumulator
    {
        public:
        
            typedef std::string Key;
        
            struct Deviation
            {
                // First moment
                double mean;
            
                // Standard deviation
                double sd;
            
                inline std::string operator()() const
                {
                    return (boost::format("%1% \u00B1 %2%") % mean % sd).str();
                }
            };

            void add(const Key &key, double value)
            {
                _data[key].push_back(value);
            }

            void add(const Key &key, const Limit &l)
            {
                _limits[key].push_back(l);
            }
        
            Deviation value(const Key &key) const
            {
                Deviation d;

                d.sd   = SS::sd(_data.at(key));
                d.mean = SS::mean(_data.at(key));
            
                return d;
            }
        
            const Limit & limits(const Key &key) const
            {
                const Limit *min = nullptr;
            
                for (const auto &limit : _limits.at(key))
                {
                    if (!min || limit.abund < min->abund)
                    {
                        min = &limit;
                    }
                }
            
                return *min;
            }
        
        private:
        
            // Used for comparing limit of detection
            std::map<Key, std::vector<Limit>> _limits;
        
            // Used for regular mappings
            std::map<Key, std::vector<double>> _data;
    };
}

#endif