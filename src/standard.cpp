#include <assert.h>
#include <algorithm>
#include "Standard.hpp"

bool Standard::known(const GeneID &id) const
{
    return std::find_if(genes.begin(), genes.end(), [&](const Gene &g)
                 {
                     return (g.id == id);
                 }) != genes.end();
}

/*
Concentration Standard::concent(const GeneID &id) const
{
    const auto iter = std::find_if(mixA.begin(), mixA.end(), [&](const MixturePair &p)
                 {
                     return (p.x.id == id || p.y.id == id);
                 });
    assert(iter != mixA.end());

    return (iter->x.id == id) ? iter->x : iter->y;
}
*/
