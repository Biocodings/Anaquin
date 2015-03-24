#ifndef AS_BUILDER_GTF_HPP
#define AS_BUILDER_GTF_HPP

#include "Feature.hpp"

/*
 * ParserGTF is simply a parser, it is unable to construct individual features to transcripts.
 * This class is designed to fit the gap and therefore complement with ParserGTF.
 */

struct BuilderGTF
{
    void add(const Feature &);
};

#endif