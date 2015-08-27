#include <iostream>
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

struct MixtureData
{
    MixtureData(const SequinID &id, Base length, Concentration abund) : id(id), length(length), abund(abund) {}
    
    inline bool operator<(const MixtureData &x)  const { return id < x.id;  }
    inline bool operator==(const MixtureData &x) const { return id == x.id; }
    
    SequinID id;
    
    // Length of the sequin
    Base length;
    
    // Amount of spiked-in concentration
    Concentration abund;
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
    _mixes  = std::shared_ptr<Mixtures>(new Mixtures);
    _annots = std::shared_ptr<Annotations>(new Annotations);
}

std::size_t Reference::countMixes() const
{
    return _mixes->size();
}

const SequinData * Reference::seq(const SequinID &id) const
{
    return _data.count(id) ? &_data.at(id) : nullptr;
}

void Reference::validate()
{
    std::vector<SequinID> mixIDs, antIDs;
    
    antIDs.resize((*_annots).size());
    mixIDs.resize((*_mixes)[MixA].size());

    std::transform((*_mixes)[MixA].begin(), (*_mixes)[MixA].end(), mixIDs.begin(), [&](const MixtureData &m)
    {
        return m.id;
    });
    
    std::transform((*_annots).begin(), (*_annots).end(), antIDs.begin(), [&](const AnnotationData &m)
    {
        return m.id;
    });

    /*
     * Check for any sequin defined in mixture but not in annotation
     */

    std::vector<SequinID> diffs, inters;

    std::set_difference(mixIDs.begin(),
                        mixIDs.end(),
                        antIDs.begin(),
                        antIDs.end(),
                        std::back_inserter(diffs));

    std::set_intersection(mixIDs.begin(),
                          mixIDs.end(),
                          antIDs.begin(),
                          antIDs.end(),
                          std::back_inserter(inters));

    /*
     * Construct a set of validated sequins
     */
    
    std::for_each(mixIDs.begin(), mixIDs.end(), [&](const SequinID &id)
    {
        auto d = SequinData();

        // The rest of the fields will be filled later...
        d.id = id;

        // Add a new entry for the validated sequin
        _data[id] = d;
    });

    /*
     * Now, we have a list of validated sequins. Use those sequins to combine information
     * from mixtures and annotations.
     */

    for (const auto i : (*_mixes))
    {
        // Eg: MixA, MixB etc
        const auto mix = i.first;

        // For each of the mixture defined
        for (const auto j : i.second)
        {
            // Only if it's a validated sequin
            if (_data.count(j.id))
            {
                _data.at(j.id).length = j.length;
                _data.at(j.id).mixes[mix] = j.abund;
            }
        }
    }
}

void Reference::add(const SequinID &id, Base length, Concentration c, Mixture m)
{
    (*_mixes)[m].insert(MixtureData(id, length, c));
}

void Reference::add(const SequinID &id, const Locus &l)
{
    _annots->insert(AnnotationData(id, l));
}