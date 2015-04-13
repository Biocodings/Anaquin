#ifndef GI_WRITER_HPP
#define GI_WRITER_HPP

#include <string>

namespace Spike
{
    struct Writer
    {
        virtual void write(const std::string &line) = 0;
    };
}

#endif