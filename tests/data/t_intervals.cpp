#include <catch.hpp>
#include "data/intervals.hpp"

#include "data/itree.hpp"

using namespace Anaquin;

TEST_CASE("Interval_Test_1")
{
    /*
     * (0,2) -> 0
     * (3,5) -> 4
     * (6,6) -> 0
     * (7,7) -> 2
     * (8,8) -> 1
     * (9,9) -> 0
     *
     *  0,0,0,0,0,1,2,4,4,4
     */
    
    Interval i("Test", Locus(0, 9));

    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(7, 7));
    i.add(Locus(7, 7));
    i.add(Locus(8, 8));

    const auto r = i.stats();
    
    REQUIRE(r.min      == 0);
    REQUIRE(r.max      == 4);
    REQUIRE(r.length   == 10);
    REQUIRE(r.zeros    == 5);
    REQUIRE(r.nonZeros == 5);
    REQUIRE(r.sums     == 15);
    REQUIRE(r.p25      == 0);
    REQUIRE(r.p50      == 0);
    REQUIRE(r.p75      == 2);
}

TEST_CASE("Interval_Test_2")
{
    /*
     * (0,4) -> 5
     * (5,5) -> 1
     * (6,9) -> 5
     *
     *  1,5,5,5,5,5,5,5,5
     */
    
    Interval i("Test", Locus(0, 9));
    
    i.add(Locus(0, 4));
    i.add(Locus(0, 4));
    i.add(Locus(0, 4));
    i.add(Locus(0, 4));
    i.add(Locus(0, 4));
    i.add(Locus(5, 5));
    i.add(Locus(6, 9));
    i.add(Locus(6, 9));
    i.add(Locus(6, 9));
    i.add(Locus(6, 9));
    i.add(Locus(6, 9));

    const auto r = i.stats();
    
    REQUIRE(r.min      == 1);
    REQUIRE(r.max      == 5);
    REQUIRE(r.length   == 10);
    REQUIRE(r.zeros    == 0);
    REQUIRE(r.nonZeros == 10);
    REQUIRE(r.sums     == 46);
    REQUIRE(r.p25      == 5);
    REQUIRE(r.p50      == 5);
    REQUIRE(r.p75      == 5);
}

TEST_CASE("Interval_Test_3")
{
    Interval i("Test", Locus(0, 9));
    
    i.add(Locus(4, 6));
    i.add(Locus(4, 6));
    i.add(Locus(4, 6));
    i.add(Locus(1, 2));
    i.add(Locus(3, 3));
    
    std::vector<Base> x, y, z;
    
    i.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
    {
        x.push_back(i);
        y.push_back(j);
        z.push_back(depth);
    });

    REQUIRE(x.size() == 4);
    REQUIRE(y.size() == 4);
    REQUIRE(z.size() == 4);
    
    REQUIRE(x[0] == 0);
    REQUIRE(x[1] == 1);
    REQUIRE(x[2] == 4);
    REQUIRE(x[3] == 7);

    REQUIRE(y[0] == 1);
    REQUIRE(y[1] == 4);
    REQUIRE(y[2] == 7);
    REQUIRE(y[3] == 10);

    REQUIRE(z[0] == 0);
    REQUIRE(z[1] == 1);
    REQUIRE(z[2] == 3);
    REQUIRE(z[3] == 0);
}






namespace Anaquin
{
    class Interval__
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
        
        Interval__(const IntervalID &id, const Locus &l) : _id(id), _l(l)
        {
            _covs.resize(l.length());
        }
        
        inline void add(const Locus &l)
        {
            if (l.start < _covs.size()) { _covs[l.start].starts++; }
            if (l.end < _covs.size())   { _covs[l.end].ends++;     }
            else                        { _covs.back().ends++;     }
        };
        
        inline Base map(const Locus &l, Base *lp = nullptr, Base *rp = nullptr)
        {
            const auto start = std::max(_l.start, l.start) - _l.start;
            const auto end   = std::min(_l.end,   l.end)   - _l.start;
            
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
            
            // Bases to the left of the interval fails to map
            const auto left  = ((l.start < _l.start) ? _l.start -  l.start : 0);
            
            // Bases to the right of the interval fails to map
            const auto right = ((l.end   > _l.end)   ?  l.end   - _l.end   : 0);
            
            if (lp) { *lp = left;  }
            if (rp) { *rp = right; }
            
            return left + right;
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
        
        // The represented interval
        Locus _l;
        
        IntervalID _id;
        
        // For each base in the interval (relative to the beginning of the interval)
        std::vector<Depth> _covs;
    };
    
    template <typename T = Interval__> class Intervals_
    {
    public:
        
        typedef std::map<typename T::IntervalID, T> IntervalData;
        
        inline void add(const T &i)
        {
            _inters.insert(typename std::map<typename T::IntervalID, T>::value_type(i.id(), i));
        }
        
        inline void build()
        {
            std::vector<Interval_<Interval__ *>> loci;
            
            #define LOCUS_TO_TINTERVAL(x) Interval_<Interval__ *>(x.l().start, x.l().end, &x)
            
            for (auto &i : _inters)
            {
                loci.push_back(LOCUS_TO_TINTERVAL(i.second));
            }
            
            assert(!loci.empty());
            
            _tree = std::shared_ptr<IntervalTree<Interval__ *>>(new IntervalTree<Interval__ *> { loci });
        }
        
        inline T * find(const typename T::IntervalID &id)
        {
            return _inters.count(id) ? &(_inters.at(id)) : nullptr;
        }
        
        inline const T * find(const typename T::IntervalID &id) const
        {
            return _inters.count(id) ? &(_inters.at(id)) : nullptr;
        }
        
        inline T * contains(const Locus &l, std::vector<T *> *r = nullptr)
        {
            auto v = _tree->findContains(l.start, l.end);

            if (r)
            {
                for (const auto &i : v)
                {
                    r->push_back(i.value);
                }
            }

            return v.empty() ? nullptr : v.front().value;
        }
        
        inline T * overlap(const Locus &l, std::vector<T *> *r = nullptr)
        {
            auto v = _tree->findOverlapping(l.start, l.end);
            
            if (r)
            {
                for (const auto &i : v)
                {
                    r->push_back(i.value);
                }
            }
            
            return v.empty() ? nullptr : v.front().value;
        }
        
        typename T::Stats stats() const
        {
            Interval::Stats stats;
            
            for (const auto &i : _inters)
            {
                const auto s = i.second.stats();
                
                stats.sums     += s.sums;
                stats.length   += s.length;
                stats.nonZeros += s.nonZeros;
                stats.zeros    += s.zeros;
                stats.min       = std::min(stats.min, s.min);
                stats.max       = std::max(stats.max, s.max);
                
                for (const auto &j : s.hist)
                {
                    stats.hist[j.first] += j.second;
                }
            }
            
            stats.mean = stats.sums / stats.length;
            
            return stats;
        }
        
        inline const IntervalData &data() const { return _inters; }
        
        // Number of intervals
        inline Counts size() const { return _inters.size(); }
        
        inline Base length() const
        {
            return std::accumulate(_inters.begin(), _inters.end(), 0,
                                   [&](int sums, const std::pair<Interval::IntervalID, Interval> & p)
                                   {
                                       return sums + p.second.l().length();
                                   });
        }
        
    private:
        
        std::shared_ptr<IntervalTree<Interval__ *>> _tree;
        
        IntervalData _inters;
    };
}

TEST_CASE("Interval_Test_4")
{
    Intervals_<> i;
    
    i.add(Interval__("1", Locus(0,  24)));
    i.add(Interval__("2", Locus(25, 49)));
    i.add(Interval__("3", Locus(50, 74)));
    i.add(Interval__("4", Locus(75, 99)));
    i.build();
    
    REQUIRE(i.contains(Locus(10, 20)));
    REQUIRE(i.contains(Locus(25, 26)));
    REQUIRE(i.contains(Locus(28, 28)));
    REQUIRE(i.contains(Locus(90, 95)));

    REQUIRE(!i.contains(Locus(00, 99)));
    REQUIRE(!i.contains(Locus(25, 50)));
    REQUIRE(!i.contains(Locus(51, 99)));
    REQUIRE(!i.contains(Locus(10, 60)));

    REQUIRE(i.overlap(Locus(00, 99)));
    REQUIRE(i.overlap(Locus(25, 50)));
    REQUIRE(i.overlap(Locus(51, 99)));
    REQUIRE(i.overlap(Locus(10, 60)));
}




