#ifndef PARSER_VCF2_HPP
#define PARSER_VCF2_HPP

#include "data/reader.hpp"
#include "data/variant.hpp"

namespace Anaquin
{
    struct ParserVCF2
    {
        typedef std::function<void (Variant &)> Functor;
        static void parse(const Reader &r, Functor f);
    };
}

#endif
