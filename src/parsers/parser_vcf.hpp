#ifndef PARSER_VCF_HPP
#define PARSER_VCF_HPP

#include "data/reader.hpp"
#include "data/variant.hpp"

namespace Anaquin
{
    struct ParserVCF
    {
        typedef std::function<void (Variant &)> Functor;
        static void parse(const Reader &r, Functor f);
    };
}

#endif
