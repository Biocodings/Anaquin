#ifndef GI_PARSER_HPP
#define GI_PARSER_HPP

#include <memory>
#include "writers/writer.hpp"
#include "exceptions/invalid_extension.hpp"

namespace Spike
{
    struct ParserOptions
    {
        std::shared_ptr<Writer> writer;
    };
}

#endif