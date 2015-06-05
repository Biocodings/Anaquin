#ifndef GI_PARSER_HPP
#define GI_PARSER_HPP

namespace Spike
{
    struct ParserProgress
    {
        long long i = 0;
    };

    enum DataMode
    {
        File,
        String,
    };
}

#endif