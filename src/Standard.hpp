#ifndef AS_STANDARD_HPP
#define AS_STANDARD_HPP

#include <set>
#include <map>
#include <list>
#include <vector>
#include "Locus.hpp"
#include "Feature.hpp"
//#include "Sequence.hpp"
#include "ConfusionMatrix.hpp"

enum Group
{
    A,
    B,
    C,
    D
};

struct Gene
{
    GeneID id;

    // Location of the gene relative to the chromosome
    Locus l;

    std::vector<Feature> exons;
    std::vector<Feature> introns;
};

struct IMixture
{
    IsoformID id;

    // Fold ratio relative to the pair
    Fold fold;
    
    // Level of expression added to the sample for this mixture
    Expression exp;
};

struct GMixture
{
    Group gr;
    GeneID id;
    
    // Reference mixture
    IMixture r;
    
    // Variant mixture
    IMixture v;
    
    // Total level of expression
    Expression exp;
};

struct Standard
{
	inline bool matchFeature(const Feature &q) const
	{
		for (auto r : fs)
		{
			if (r.type == q.type && r.l.contains(q.l))
			{
				return true;
			}
		}

		return false;
	}
    
    // Whether the given gene is a part of the standard
    bool known(const GeneID &id) const;
    
    ChromoID id;

    // The location of the chromosome
    Locus l;

    std::map<GeneID, GMixture> mixA;
    std::map<GeneID, GMixture> mixB;

    std::map<IsoformID, IMixture> isoA;
    std::map<IsoformID, IMixture> isoB;

    std::map<IsoformID, GeneID> iso2Gene;
    
    // Known genes
    std::vector<Gene> genes;

    // Known features
    std::vector<Feature> fs;

    // Known exons
    std::vector<Feature> exons;

    // Known introns (spliced junctions)
    std::vector<Feature> introns;
};

#endif