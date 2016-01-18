#ifndef MOCK_WRITER_HPP
#define MOCK_WRITER_HPP

#include "writers/writer.hpp"

namespace Anaquin
{
    struct MockWriter : public Writer
    {
        inline void close() override {}
        inline void open(const FileName &) override {}
        inline void create(const std::string &) override {}
        inline void write(const std::string &, bool newLine) override {}
    };
}

#endif