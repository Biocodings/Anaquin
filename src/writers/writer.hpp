#ifndef WRITER_HPP
#define WRITER_HPP

#include <string>

namespace Anaquin
{
    template <typename T = std::string> struct Writer
    {
        virtual void close() = 0;
        virtual void open(const FileName &) = 0;
        virtual void write(const T &) = 0;
    };
}

#endif
