#ifndef WRITER_HPP
#define WRITER_HPP

#include <string>

namespace Anaquin
{
    struct Writer
    {
        virtual void close() = 0;
        virtual void open(const FileName &file) = 0;
        virtual void write(const std::string &line) = 0;
    };
}

#endif