#include <iostream>
#include <assert.h>
#include "m_align.hpp"
#include "biology.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "writers/writer.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

MAlignStats MAlign::analyze(const std::string &file, const MAlign &options)
{
	return MAlignStats();
}