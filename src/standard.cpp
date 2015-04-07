#include <assert.h>
#include <algorithm>
#include "standard.hpp"

using namespace Spike;

bool Standard::known(const GeneID &id) const
{
    return std::find_if(genes.begin(), genes.end(), [&](const Gene &g)
                 {
                     return (g.id == id);
                 }) != genes.end();
}