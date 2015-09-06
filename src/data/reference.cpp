#include "data/tokens.hpp"
#include "data/reference.hpp"

#include <iostream>

using namespace Anaquin;

struct VarRef::VariantPair
{
    const MixtureData *r, *v;
};

struct VarRef::VarRefImpl
{
    /*
     * Validated variables
     */

    std::vector<Variation> vars;

    std::map<Mixture, std::map<BaseID, VariantPair>> pairs;
    
    /*
     * Raw variables
     */
    
    std::vector<Variation> rawVars;
};

VarRef::VarRef() : _impl(new VarRefImpl()) {}

std::size_t VarRef::countVars() const
{
    return _impl->vars.size();
}

double VarRef::alleleFreq(Mixture m, const BaseID &bID) const
{
    const auto &p = _impl->pairs.at(m).at(bID);
    const auto &r = p.r;
    const auto &v = p.v;

    // Abundance ratio of reference to variant DNA standard
    return v->abund / (r->abund + v->abund);
}

void VarRef::addVar(const Variation &v)
{
    assert(!v.bID.empty());
    _impl->rawVars.push_back(v);
}

void VarRef::validate()
{
    _impl->vars = _impl->rawVars;

    // Validate sequins defined in the mixture
    merge(_rawMIDs, _rawMIDs);
    
   /*
    * Construct data structure for homozygous/heterozygous
    */

    std::vector<std::string> toks;
    
    for (const auto &i : _mixes)
    {
        for (const auto &j : _mixes.at(i.first))
        {
            Tokens::split(j.id, "_", toks);

            // It has be the reference or variant...
            assert(toks[3] == "R" || toks[3] == "V");
            
            // Eg: D_1_10
            const auto baseID = toks[0] + "_" + toks[1] + "_" + toks[2];

            if (toks[3] == "R")
            {
                _impl->pairs[i.first][baseID].r = &j;
            }
            else
            {
                _impl->pairs[i.first][baseID].v = &j;
            }
        }
    }

    assert(!_impl->pairs.empty());
}

const Variation * VarRef::findVar(const Locus &l, double fuzzy, Matching match) const
{
    for (const auto &i : _impl->vars)
    {
        assert(i.l.length());

        if (match == StartOnly && i.l.start == l.start)
        {
            return &i;
        }
    }

    return nullptr;
}