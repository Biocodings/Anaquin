#ifndef INTERVALS_HPP
#define INTERVALS_HPP

#include <map>
#include <numeric>
#include "data/locus.hpp"
#include "data/types.hpp"

namespace Anaquin
{
    class Interval
    {
        public:
        
            typedef std::string IntervalID;
        
            struct Stats
            {
                Coverage min = std::numeric_limits<Coverage>::max();
                Coverage max = std::numeric_limits<Coverage>::min();

                // Percentile for the interval
                Coverage p25, p50, p75;
                
                // Arithmetic first moment
                double mean;
                
                // Sum of all coverage
                Coverage sums = 0;

                // Distribution for the coverage
                std::map<Coverage, Counts> hist;
                
                // Length of the interval
                Base length = 0;

                // Number of bases with non-zero coverage
                Base nonZeros = 0;

                // Number of bases with zero coverage
                Counts zeros = 0;

                inline double covered() const { return static_cast<double>(nonZeros) / length; }
            };

            Interval(const IntervalID &id, const Locus &l) : _id(id), _l(l)
            {
                _covs.resize(l.length());
            }

            inline void add(const Locus &l)
            {
                if (l.start < _covs.size()) { _covs[l.start].starts++; }
                if (l.end < _covs.size())   { _covs[l.end].ends++;     }
                else                        { _covs.back().ends++;     }
            };

            inline Base add_(const Locus &l)
            {
                auto start = std::max(_l.start, l.start) - _l.start;
                auto end   = std::min(_l.end,   l.end)   - _l.start;

                /*
                 * For example, if the interval is (2,2038) and the locus is (2028,2042).
                 *
                 *     start = 2028 -> 2026
                 *     end   = 2038 -> 2036
                 */
                
                if (start <= end)
                {
                    _covs[start].starts++;
                    _covs[end].ends++;
                    
                    assert(start < _covs.size() && end < _covs.size());
                }
                
                return ((l.start < _l.start) ? _l.start -  l.start : 0) +
                       ((l.end   > _l.end)   ?  l.end   - _l.end   : 0);
            }

            template <typename T> Stats stats(T t) const
            {
                Stats stats;

                bedGraph([&](const ChromoID &id, Base i, Base j, Coverage cov)
                {
                    // Should this be counted? For example, aligning to sequins?
                    if (!t(id, i, j, cov))
                    {
                        return;
                    }

                    stats.min = std::min(stats.min, cov);
                    stats.max = std::max(stats.max, cov);
                    
                    // The interval is half-open
                    const auto n = j - i;

                    stats.sums      += (n * cov);
                    stats.length    += n;
                    stats.hist[cov] += n;

                    if (!cov) { stats.zeros    += n; }
                    else      { stats.nonZeros += n; }
                });

                auto percent = [&](Counts n)
                {
                    Counts i = 0;
                    
                    for (const auto &p : stats.hist)
                    {
                        i += p.second;
                        
                        // Have we reached our limit?
                        if (i >= n)
                        {
                            return p.first;
                        }
                    }

                    return stats.max;
                };

                stats.mean = stats.sums / stats.length;
                stats.p25  = percent(0.25 * stats.length);
                stats.p50  = percent(0.50 * stats.length);
                stats.p75  = percent(0.75 * stats.length);

                return stats;
            }

            inline Stats stats() const
            {
                return stats([&](const ChromoID &id, Base i, Base j, Coverage cov)
                {
                    return true;
                });
            }
        
            template <typename T> void bedGraph(T t) const
            {
                Base depth = 0;
                long lastStart = -1;
                long lastDepth = -1;

                for (auto j = 0; j < _l.length(); j++)
                {
                    depth += _covs[j].starts;
                
                    if (depth != lastDepth)
                    {
                        if (lastDepth != -1)
                        {
                            t(_id, lastStart, j, lastDepth);
                        }

                        // Set current position as the new interval start + depth
                        lastStart = j;
                        lastDepth = depth;
                    }

                    depth = depth - _covs[j].ends;
                }
            
                // Print information about the last position
                if (lastDepth != -1)
                {
                    t(_id, lastStart, _l.length(), lastDepth);
                }
            }

            inline Locus l()       const { return _l;  }
            inline IntervalID id() const { return _id; }

        private:

            struct Depth
            {
                Base starts;
                Base ends;
            };
        
            Locus _l;
            IntervalID _id;

            // For each base in the interval
            std::vector<Depth> _covs;
    };

    class Intervals
    {
        public:

            typedef std::map<Interval::IntervalID, Interval> IntervalMap;

            inline void add(const Interval &i)
            {
                _inters.insert(std::map<Interval::IntervalID, Interval>::value_type(i.id(), i));
            }

            inline Interval * find(const Interval::IntervalID &id)
            {
                return _inters.count(id) ? &(_inters.at(id)) : nullptr;
            }

            inline const Interval * find(const Interval::IntervalID &id) const
            {
                return _inters.count(id) ? &(_inters.at(id)) : nullptr;
            }

            inline Interval * contains(const Locus &l)
            {
                for (auto &i : _inters)
                {
                    if (i.second.l().contains(l))
                    {
                        return &i.second;
                    }
                }

                return nullptr;
            }

            inline Interval *overlap(const Locus &l)
            {
                for (auto &i : _inters)
                {
                    if (i.second.l().overlap(l))
                    {
                        return &i.second;
                    }
                }

                return nullptr;
            }
        
            Interval::Stats stats() const;
        
            inline const IntervalMap &map() const { return _inters; }

            // Number of intervals
            inline std::size_t size() const { return _inters.size(); }

            inline Base length() const
            {
                return std::accumulate(_inters.begin(), _inters.end(), 0,
                          [&](int sums, const std::pair<Interval::IntervalID, Interval> & p)
                          {
                              return sums + p.second.l().length();
                          });
            }

        private:
            IntervalMap _inters;
    };
}

#endif