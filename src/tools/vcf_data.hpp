#ifndef VCF_DATA_HPP
#define VCF_DATA_HPP

#include "data/hist.hpp"
#include "data/standard.hpp"
#include "data/intervals.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser_vcf.hpp"

namespace Anaquin
{
    typedef long VarKey;
    
    struct VCFChrData
    {
        // SNP to data
        std::map<Base, Variant> s2d;
        
        // Indels to data
        std::map<Base, Variant> i2d;
    };

    struct VCFData : public std::map<ChrID, VCFChrData>
    {
        inline std::map<ChrID, std::map<long, Counts>> hist() const
        {
            std::map<ChrID, std::map<long, Counts>> r;
            
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.s2d)
                {
                    const auto key = var2hash(j.second.id, j.second.type(), j.second.l);
                    r[i.first][key];
                }

                for (const auto &j : i.second.i2d)
                {
                    const auto key = var2hash(j.second.id, j.second.type(), j.second.l);
                    r[i.first][key];
                }
                
                assert(!r[i.first].empty());
            }
            
            assert(!r.empty());
            return r;
        }

        inline const Variant * findVar(const ChrID &cID, VarKey key)
        {
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.s2d)
                {
                    if (key == var2hash(j.second.id, j.second.type(), j.second.l))
                    {
                        return &j.second;
                    }
                }
                
                for (const auto &j : i.second.i2d)
                {
                    if (key == var2hash(j.second.id, j.second.type(), j.second.l))
                    {
                        return &j.second;
                    }
                }
            }

            return nullptr;
        }
        
        inline const Variant * findVar(const ChrID &cID, const Locus &l)
        {
            if (at(cID).s2d.count(l.start))
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
            return at(cID).s2d.size();
        }
        
        inline Counts countVar() const
        {
            return countSNP() + countInd();
        }
        
        inline Counts countSNP() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const VCFChrData &x)
            {
                return countSNP(cID);
            });
        }
        
        inline Counts countSNPSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const VCFChrData &x)
            {
                return Standard::isSynthetic(cID) ? countSNP(cID) : 0;
            });
        }
        
        inline Counts countSNPGen() const
        {
            return countSNP() - countSNPSyn();
        }
        
        inline Counts countInd(const ChrID &cID) const
        {
            return at(cID).i2d.size();
        }
        
        inline Counts countInd() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const VCFChrData &x)
            {
                return countInd(cID);
            });
        }
        
        inline Counts countIndSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const VCFChrData &x)
            {
                return Standard::isSynthetic(cID) ? countInd(cID) : 0;
            });
        }
    
        inline Counts countIndGen() const
        {
            return countInd() - countIndSyn();
        }
    };
    
    inline VCFData vcfData(const Reader &r)
    {
        VCFData c2d;
        
        ParserVCF::parse(r, [&](const ParserVCF::Data &x, const ParserProgress &)
        {
            switch (x.type())
            {
                case Mutation::SNP:
                {
                    c2d[x.cID].s2d[x.l.start] = x;
                    break;
                }

                case Mutation::Deletion:
                case Mutation::Insertion:
                {
                    c2d[x.cID].i2d[x.l.start] = x;
                    break;
                }
            }
        });

        return c2d;
    }    
}

#endif
