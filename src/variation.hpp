#ifndef GI_VARIATION_HPP
#define GI_VARIATION_HPP

#include "types.hpp"
#include "sequins.hpp"

struct VariationStats
{
};

struct Variation
{
	static VariationStats analyze(const std::string &file);
};

#endif