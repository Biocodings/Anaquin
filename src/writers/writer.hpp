#ifndef GI_WRITER_HPP
#define GI_WRITER_HPP

#include <string>

namespace Anaquin
{
    struct Writer
    {
        virtual void close() = 0;
        virtual void open(const std::string &file) = 0;
        virtual void write(const std::string &line) = 0;
    };
}

#endif