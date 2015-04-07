#ifndef GI_READER_FA_HPP
#define GI_READER_FA_HPP

#include <functional>
#include "types.hpp"

namespace Spike
{
    struct FASequence
    {
        std::string id;
        Sequence value;
    };
    
    struct ParserFA
    {
        static void parse(const std::string &file, std::function<void (const FASequence &)>);
    };    
}

#endif