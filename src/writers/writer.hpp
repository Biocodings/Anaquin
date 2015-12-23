#ifndef WRITER_HPP
#define WRITER_HPP

#include <string>

namespace Anaquin
{
    struct Writer
    {
        virtual void close() = 0;
        virtual void open(const FileName &) = 0;
        virtual void write(const std::string &) = 0;
        virtual void create(const std::string &) = 0;
    };
}

#endif