#ifndef MERGED_HPP
#define MERGED_HPP

#include <map>
#include <numeric>
#include "data/data.hpp"
#include "data/itree.hpp"
#include "data/locus.hpp"

namespace Anaquin
{
    class MergedInterval : public Matched
    {
        public:
        
            typedef std::string IntervalID;
        
            struct Stats
            {
                Base length = 0;
                Base nonZeros = 0;
            
                inline Proportion covered() const { return static_cast<double>(nonZeros) / length; }
            };
        
            MergedInterval(const IntervalID &id, const Locus &l) : _id(id), _l(l)
            {
                _id = id;
            }

            inline Base map(const Locus &l, Base *lp = nullptr, Base *rp = nullptr)
            {
                bool added = true;
                bool p1 = true;
                bool p2 = false;
                
                // Pointing to the matching
                Locus *m = nullptr;

                Base left = 0;
                Base right = 0;
                
                for (auto it = _data.begin(); it != _data.cend();)
                {
                    auto &j = it->second;
                    
                    if (p1)
                    {
                        if (j.overlap(l))
                        {
                            added = false;
                            
                            left  = ((l.start < j.start) ? j.start -  l.start : 0);
                            right = ((l.end   > j.end)   ?  l.end  - j.end   : 0);
                            
                            if (right > 10)
                            {
                                right = right;
                            }
                            

                            j.start = std::min(j.start, l.start);
                            j.end   = std::max(j.end,   l.end);
                            
                            j.start = std::max(j.start, _l.start);
                            j.end   = std::min(j.end,   _l.end);
                            
                            
                            // So that we can access it in the later iterations
                            m = &j;

                            if (lp) { *lp = left;  }
                            if (rp) { *rp = right; }
                            
                            p1 = false;
                            p2 = true;
                        }
                    }
                    else if (p2)
                    {
                        if (!j.overlap(l))
                        {
                            return (left + right);
                        }
                        
                        m->start = std::min(m->start, j.start);
                        m->end   = std::max(m->end,   j.end);
                        
                        m->start = std::max(m->start, _l.start);
                        m->end   = std::min(m->end,   _l.end);
                        
                        // Remove the overlapping entry
                        _data.erase(it++);
                        
                        continue;
                    }
                    else
                    {
                        throw "Invalid phase";
                    }
                    
                    ++it;
                }
                
                if (added)
                {
                    auto start = l.start;
                    auto end = l.end;
                    
                    start = std::max(start, _l.start);
                    end   = std::min(end,   _l.end);
                    
                    _data[start] = Locus(start, end);
                    
                    left  = ((l.start < _l.start) ? _l.start -  l.start : 0);
                    right = ((l.end   > _l.end)   ?  l.end  - _l.end   : 0);
                }
                
                if (lp) { *lp = left;  }
                if (rp) { *rp = right; }
                
                return left + right;
            }
        
            template <typename F> Stats stats(F f) const
            {
                Stats stats;

                for (const auto &i : _data)
                {
                    stats.nonZeros += i.second.length();
                }
                
                stats.length = _l.length();
                assert(stats.length >= stats.nonZeros);

                return stats;
            }
        
            inline Stats stats() const
            {
                return stats([&](const ChrID &id, Base i, Base j, Coverage cov)
                {
                    return true;
                });
            }
        
            inline const Locus &l()       const { return _l;  }
            inline const IntervalID &id() const { return _id; }
        
            inline IntervalID name() const override { return id(); }
        
            inline std::size_t size() { return _data.size(); }
        
        private:
        
            Locus _l;

            std::map<Base, Locus> _data;
        
            IntervalID _id;
    };
    
    template <typename T = MergedInterval> class MergedIntervals
    {
        public:
        
            typedef std::map<typename T::IntervalID, T> IntervalData;

            inline void add(const T &i)
            {
                _inters.insert(typename std::map<typename T::IntervalID, T>::value_type(i.id(), i));
            }

            inline void build()
            {
                std::vector<Interval_<T *>> loci;
            
                #define LOCUS_TO_TINTERVAL(x) Interval_<T *>(x.l().start, x.l().end, &x)
            
                for (auto &i : _inters)
                {
                    loci.push_back(LOCUS_TO_TINTERVAL(i.second));
                }
                
                if (loci.empty())
                {
                    throw "No interval was built. loci.empty().";
                }
            
                _tree = std::shared_ptr<IntervalTree<T *>>(new IntervalTree<T *> { loci });
            }
        
            inline T * find(const typename T::IntervalID &id)
            {
                return _inters.count(id) ? &(_inters.at(id)) : nullptr;
            }
        
            inline const T * find(const typename T::IntervalID &id) const
            {
                return _inters.count(id) ? &(_inters.at(id)) : nullptr;
            }

            inline T * exact(const Locus &l, std::vector<T *> *r = nullptr) const
            {
                auto v = _tree->findContains(l.start, l.end);

                T *t = nullptr;
    
                for (const auto &i : v)
                {
                    if (i.value->l() == l)
                    {
                        t = i.value;
                        
                        if (r)
                        {
                            r->push_back(t);
                        }
                    }
                }
            
                return t;
            }
        
            inline T * contains(const Locus &l, std::vector<T *> *r = nullptr) const
            {
                auto v = _tree->findContains(l.start, l.end);

                if (r)
                {
                    for (const auto &i : v)
                    {
                        if (i.value)
                        
                        r->push_back(i.value);
                    }
                }
            
                return v.empty() ? nullptr : v.front().value;
            }
        
            inline T * overlap(const Locus &l, std::vector<T *> *r = nullptr) const
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
                MergedInterval::Stats stats;
            
                for (const auto &i : _inters)
                {
                    const auto s = i.second.stats();
                
                    stats.length   += s.length;
                    stats.nonZeros += s.nonZeros;
                }
            
                return stats;
            }
        
            inline const IntervalData &data() const { return _inters; }
        
            // Number of intervals
            inline Counts size() const { return _inters.size(); }
        
            inline Base length() const
            {
                return std::accumulate(_inters.begin(), _inters.end(), 0,
                        [&](int sums, const std::pair<MergedInterval::IntervalID, MergedInterval> & p)
                {
                    return sums + p.second.l().length();
                });
            }
        
        private:
        
            std::shared_ptr<IntervalTree<T *>> _tree;
        
            IntervalData _inters;
    };

    typedef std::map<ChrID, MergedIntervals<>> MC2Intervals;

    struct MergedID2Intervals : std::map<MergedInterval::IntervalID, MergedIntervals<>>
    {
        inline void add(const MergedInterval::IntervalID &id, const MergedIntervals<> &i)
        {
            (*this)[id] = i;
        }
        
        inline Counts countInters() const
        {
            Counts n = 0;
            
            for (const auto &i : *this)
            {
                n += i.second.size();
            }
            
            return n;
        }
        
        inline MergedInterval::Stats stats() const
        {
            MergedInterval::Stats stats;
            
            for (const auto &i : *this)
            {
                const auto s = i.second.stats();
                
                stats.length   += s.length;
                stats.nonZeros += s.nonZeros;
            }
            
            return stats;
        }
    };
}

#endif