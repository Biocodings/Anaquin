#ifndef VCF_DATA_HPP
#define VCF_DATA_HPP

#include "tools/tools.hpp"
//#include "data/standard.hpp"
#include "parsers/parser_vcf.hpp"

namespace Anaquin
{
    typedef long VarKey;
    
    struct VDataData
    {
        std::map<Base, Variant> b2v;
        std::map<Variation, std::set<Variant>> m2v;
    };

    struct VCFData : public std::map<ChrID, VDataData>
    {
        inline std::set<Variant> vars() const
        {
            std::set<Variant> x;
            
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.m2v)
                {
                    for (const auto &k : j.second)
                    {
                        x.insert(k);
                    }
                }
            }
            
            return x;
        }

        inline const Variant * findVar(const ChrID &id, const Locus &l)
        {
            if (!count(id))
            {
                return nullptr;
            }
            else if (at(id).b2v.count(l.start))
            {
                return &(at(id).b2v.at(l.start));
            }

            return nullptr;
        }

        inline Counts count_(const ChrID &cID, Variation m) const
        {
            return count(cID) && at(cID).m2v.count(m) ? at(cID).m2v.at(m).size() : 0;
        }

        inline Counts count_(Variation m) const
        {
            return countMap(*this, [&](const ChrID &cID, const VDataData &)
            {
                return count_(cID, m);
            });
        }
    };

    template <typename F> VCFData readVFile(const Reader &r, F f)
    {
        VCFData c2d;
        
        ParserVCF::parse(r, [&](const Variant &x)
        {
            c2d[x.cID].b2v[x.l.start] = x;
            c2d[x.cID].m2v[x.type()].insert(x);
            f(x);
        });
        
        return c2d;
    }
    
    inline VCFData readVFile(const Reader &r)
    {
        return readVFile(r, [](const Variant &) {});
    }
}

#endif
