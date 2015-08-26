#ifndef GI_REFERENCE_HPP
#define GI_REFERENCE_HPP

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

    struct MixtureData
    {
        MixtureData(const SequinID &id, Base length, Concentration c) : id(id), length(length), c(c) {}
        
        inline bool operator<(const MixtureData &x)  const { return id < x.id;  }
        inline bool operator==(const MixtureData &x) const { return id == x.id; }
        
        SequinID id;
        
        // Length of the sequin
        Base length;
        
        // Raw amount of spiked-in concentration
        Concentration c;
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

            // Returns the valid sequins
            const std::set<SequinID> &seqs() const;

            const MixtureData &seq(const SequinID &id, Mixture m) const;

            void valid();

        private:

            struct Stats;
            struct Mixtures;
            struct Annotations;

            std::shared_ptr<Stats>       _stats;
            std::shared_ptr<Mixtures>    _rawMixes;
            std::shared_ptr<Annotations> _rawAnnots;
            std::shared_ptr<Mixtures>    _validMixes;
            std::shared_ptr<Annotations> _validAnnots;
    };
}

#endif