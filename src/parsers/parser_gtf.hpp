#ifndef GI_PARSER_GTF_HPP
#define GI_PARSER_GTF_HPP

#include "feature.hpp"
#include "analyzer.hpp"

namespace Spike
{
    struct ParserGTF
    {
        typedef std::function<void (const Feature &, const ParserProgress &)> Callback;
        static void parse(const std::string &, Callback, DataMode mode = File);
    };
}

#endif