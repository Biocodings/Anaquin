#ifndef GI_REFERENCE_HPP
#define GI_REFERENCE_HPP

#include <set>
#include <memory>

namespace Anaquin
{
    enum Mixture
    {
        MixA,
        MixB,
        MixF,
        MixG,
    };

    struct Stats
    {
    };

    struct SequinData
    {
        inline bool operator<(const SequinID &x)  const { return this->id < x;  }
        inline bool operator==(const SequinID &x) const { return this->id == x; }

        SequinID id;

        // Length of the sequin
        Base length;
        
        // Amount of spiked-in concentration
        std::map<Mixture, Concentration> mixes;

        Locus l;
    };

    class Reference
    {
        public:
            Reference();

            // Adds a sequin defined in a mixture file
            void add(const SequinID &id, Base length, Concentration c, Mixture m);

            // Adds a sequin defined in an annotation file
            void add(const SequinID &id, const Locus &l);

            // Returns number of sequins in the mixture
            std::size_t countMixes() const;

            // Returns all validated sequins
            const std::map<SequinID, SequinData> &data() const { return _data; }

            const SequinData *seq(const SequinID &id) const;

            void validate();

        private:

            struct Mixtures;
            struct Annotations;

            // Validated sequins
            std::map<SequinID, SequinData> _data;

            // Statistics about the sequins, only valid after validate()
            Stats _stats;

            std::shared_ptr<Mixtures>    _mixes;
            std::shared_ptr<Annotations> _annots;
    };
}

#endif