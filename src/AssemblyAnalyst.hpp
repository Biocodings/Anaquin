#ifndef AS_ASSEMBLY_ANALYST_HPP
#define AS_ASSEMBLY_ANALYST_HPP

#include "Sequins.hpp"
#include "ConfusionMatrix.hpp"

struct AssemblyStats
{
	ConfusionMatrix base;
	ConfusionMatrix exon;
	ConfusionMatrix intron;
};

struct AssemblyAnalyst
{
	static AssemblyStats analyze(const std::string &file, Sequins s = Sequins(), Reads n = std::numeric_limits<Reads>::max());
};

#endif