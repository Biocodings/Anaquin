#ifndef V_DATA_HPP
#define V_DATA_HPP

#include "data/hist.hpp"
#include "tools/tools.hpp"
#include "data/standard.hpp"
#include "data/intervals.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_variants.hpp"

namespace Anaquin
{
    typedef long VarKey;
    
    /*
     * Default internal representation for a set of variants
     */

    struct DefaultVarData
    {
        // SNP to data
        std::map<Base, Variant> s2d;
        
        // Indels to data
        std::map<Base, Variant> i2d;
        
        std::map<Variation, std::set<Variant>> m2v;
    };

    struct VData : public std::map<ChrID, DefaultVarData>
    {
        inline std::set<Variant> vars() const
        {
            std::set<Variant> x;
            
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.s2d)
                {
                    x.insert(j.second);
                }
                
                for (const auto &j : i.second.i2d)
                {
                    x.insert(j.second);
                }
            }
            
            return x;
        }
        
        inline std::map<ChrID, std::map<long, Counts>> hist() const
        {
            std::map<ChrID, std::map<long, Counts>> r;
            
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.s2d)
                {
                    const auto key = var2hash(j.second.name, j.second.type(), j.second.l);
                    r[i.first][key];
                }

                for (const auto &j : i.second.i2d)
                {
                    const auto key = var2hash(j.second.name, j.second.type(), j.second.l);
                    r[i.first][key];
                }
                
                assert(!r[i.first].empty());
            }
            
            assert(!r.empty());
            return r;
        }

        inline const Variant * findVar(const ChrID &cID, VarKey key)
        {
            if (!count(cID))
            {
                return nullptr;
            }
            
            //for (const auto &i : this->at(cID))
            {
                for (const auto &j : this->at(cID).s2d)
                {
                    if (key == var2hash(j.second.name, j.second.type(), j.second.l))
                    {
                        return &j.second;
                    }
                }
                
                for (const auto &j : this->at(cID).i2d)
                {
                    if (key == var2hash(j.second.name, j.second.type(), j.second.l))
                    {
                        return &j.second;
                    }
                }
            }

            return nullptr;
        }
        
        inline const Variant * findVar(const ChrID &cID, const Locus &l)
        {
            if (!count(cID))
            {
                return nullptr;
            }
            else if (at(cID).s2d.count(l.start))
            {
                return &(at(cID).s2d.at(l.start));
            }
            else if (at(cID).i2d.count(l.start))
            {
                return &(at(cID).i2d.at(l.start));
            }

            return nullptr;
        }

        inline Counts countSNP(const ChrID &cID) const
        {
            return count(cID) ? at(cID).s2d.size() : 0;
        }
        
        inline Counts countVar() const
        {
            return countSNP() + countInd();
        }
        
        inline Counts countSNP() const
        {
            return countMap(*this, [&](const ChrID &cID, const DefaultVarData &x)
            {
                return countSNP(cID);
            });
        }
        
        inline Counts countSNPSyn() const
        {
            return countMap(*this, [&](const ChrID &cID, const DefaultVarData &x)
            {
                return countSNP(cID);
            });
        }
        
        inline Counts countSNPGen() const
        {
            return countSNP() - countSNPSyn();
        }
        
        inline Counts count_(const ChrID &cID, Variation m) const
        {
            return count(cID) && at(cID).m2v.count(m) ? at(cID).m2v.at(m).size() : 0;
        }

        inline Counts count_(Variation m) const
        {
            return countMap(*this, [&](const ChrID &cID, const DefaultVarData &)
            {
                return count_(cID, m);
            });
        }

        inline Counts countInd(const ChrID &cID) const
        {
            return count(cID) ? at(cID).i2d.size() : 0;
        }
        
        inline Counts countInd() const
        {
            return countMap(*this, [&](const ChrID &cID, const DefaultVarData &x)
            {
                return countInd(cID);
            });
        }
        
        inline Counts countIndSyn() const
        {
            return countMap(*this, [&](const ChrID &cID, const DefaultVarData &x)
            {
                return countInd(cID);
            });
        }
    
        inline Counts countIndGen() const
        {
            return countInd() - countIndSyn();
        }
        
        inline Counts countVarGen() const
        {
            return countSNPGen() + countIndGen();
        }
        
        inline Counts countVarSyn() const
        {
            return countSNPSyn() + countIndSyn();
        }
    };

    template <typename F> VData readVFile(const Reader &r, F f)
    {
        VData c2d;
        
        ParserVCF::parse(r, [&](const ParserVCF::Data &x, const ParserProgress &p)
        {
            switch (x.type())
            {
                case Variation::SNP:
                {
                    c2d[x.cID].s2d[x.l.start] = x;
                    break;
                }
                    
                case Variation::Deletion:
                case Variation::Insertion:
                {
                    c2d[x.cID].i2d[x.l.start] = x;
                    break;
                }
            }

            c2d[x.cID].m2v[x.type()].insert(x);
            f(x, p);
        });
        
        return c2d;
    }
    
    inline VData readVFile(const Reader &r)
    {
        return readVFile(r, [](const ParserVCF::Data &x, const ParserProgress &) {});
    }
}

#endif
