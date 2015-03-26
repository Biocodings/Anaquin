#ifndef GI_ASSEMBLY_ANALYST_HPP
#define GI_ASSEMBLY_ANALYST_HPP

#include "sequins.hpp"
#include "confusion_matrix.hpp"

struct AssemblyStats
{
	ConfusionMatrix base;
	ConfusionMatrix exon;
	ConfusionMatrix intron;
};

struct Assembly
{
	static AssemblyStats analyze(const std::string &file, Sequins s = Sequins(), Reads n = std::numeric_limits<Reads>::max());
};

#endif