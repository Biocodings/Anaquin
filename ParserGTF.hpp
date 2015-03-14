#ifndef AS_PARSER_GTF_HPP
#define AS_PARSER_GTF_HPP

#include "Feature.hpp"

struct FeatureReader
{
    virtual void all(const Feature &)  {}
    virtual void exon(const Feature &) {}
};

struct ParserGTF
{
	static bool parse(const std::string &file, FeatureReader &reader);
};

#endif