#include <set>
#include "data/locus.hpp"
#include "data/reference.hpp"

using namespace Anaquin;

struct AnnotationData
{
    AnnotationData(const SequinID &id, Locus l) : id(id), l(l) {}

    inline bool operator<(const AnnotationData &x)  const { return id < x.id;  }
    inline bool operator==(const AnnotationData &x) const { return id == x.id; }

    SequinID id;
    Locus l;
};

struct Reference::Stats
{
    std::set<SequinID> ids;
};

struct Reference::Mixtures : public std::map<Mixture, std::set<MixtureData>>
{
    // Empty Implementation
};

struct Reference::Annotations : std::set<AnnotationData>
{
    // Empty Implementation
};

Reference::Reference()
{
    _rawMixes  = std::shared_ptr<Mixtures>(new Mixtures);
    _rawAnnots = std::shared_ptr<Annotations>(new Annotations);
}

std::size_t Reference::countMixes() const
{
    return _rawMixes->size();
}

const MixtureData & Reference::seq(const SequinID &id, Mixture m) const
{
    return *(std::find_if((*_validMixes)[m].begin(), (*_validMixes)[m].end(), [&](const MixtureData &d)
    {
        return id == d.id;
    }));
}

void Reference::valid()
{
    _stats = std::shared_ptr<Stats>(new Stats());

    for (const auto &i : (*_rawMixes)[MixA])
    {
        _stats->ids.insert(i.id);
    }
    
    _validMixes = _rawMixes; // TODO: Fix this
    
    assert(!_stats->ids.empty());
}

const std::set<SequinID> &Reference::seqs() const
{
    return _stats->ids;
}

void Reference::add(const SequinID &id, Base length, Concentration c, Mixture m)
{
    (*_rawMixes)[m].insert(MixtureData(id, length, c));
}

void Reference::add(const SequinID &id, const Locus &l)
{
    _rawAnnots->insert(AnnotationData(id, l));
}