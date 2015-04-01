#ifndef GI_DNA_VARIATION_HPP
#define GI_DNA_VARIATION_HPP

#include "types.hpp"
#include "sequins.hpp"

struct VariationStats
{
};

struct DNAVariation
{
	static VariationStats analyze(const std::string &file);
};

#endif