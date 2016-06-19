#ifndef VCF_DATA_HPP
#define VCF_DATA_HPP

#include "data/hist.hpp"
#include "data/standard.hpp"
#include "data/intervals.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser_vcf.hpp"

namespace Anaquin
{
    struct VCFChrData
    {
        std::map<Base, Variant> s2d;
        std::map<Base, Variant> i2d;
    };

    struct VCFData : public std::map<ChrID, VCFChrData>
    {
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
