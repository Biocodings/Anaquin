#ifndef MOCK_WRITER_HPP
#define MOCK_WRITER_HPP

#include "writers/writer.hpp"

namespace Anaquin
{
    struct MockWriter : public Writer
    {
        inline void close() {}
        inline void open(const FileName &file) {}
        inline void write(const std::string &line) override {}
    };
}

#endif