#ifndef MOCK_WRITER_HPP
#define MOCK_WRITER_HPP

#include "writers/writer.hpp"

namespace Anaquin
{
    struct MockWriter : public Writer
    {
        inline void close() override {}
        inline void open(const FileName &file) override {}
        inline void write(const std::string &line) override {}
    };
}

#endif