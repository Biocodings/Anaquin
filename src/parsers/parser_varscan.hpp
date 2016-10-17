#ifndef PARSER_VARSCAN_HPP
#define PARSER_VARSCAN_HPP

#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "data/convert.hpp"
#include "data/variant.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserVarScan
    {
        typedef CalledVariant Data;

        typedef std::function<void(const Data &, const ParserProgress &)> Functor;

        static bool isVarScan(const Reader &r)
        {
            return ParserVarScan::isPileup(r) || ParserVarScan::isSomatic(r);
        }

        static bool isPileup(const Reader &);
        static bool isSomatic(const Reader &);
        
        static void parsePile(const Reader &, Functor);
        static void parseSomatic(const Reader &, Functor);
        
        static void parse(const Reader &, Functor);
    };
}

#endif
