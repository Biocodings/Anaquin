#include <iostream>
#include "m_align.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "writers/writer.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

MAlignStats MAlign::analyze(const std::string &file, const Options &options)
{
	return MAlignStats();
}