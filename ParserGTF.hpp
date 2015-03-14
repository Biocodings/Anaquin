#ifndef AS_PARSER_GTF_HPP
#define AS_PARSER_GTF_HPP

#include <functional>
#include "feature.hpp"

struct ReaderGTF
{
	virtual void exon(const Feature &) = 0;
};

struct ParserGTF
{
	static bool parse(const std::string &file, ReaderGTF &reader);
};

#endif