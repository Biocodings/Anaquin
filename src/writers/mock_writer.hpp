#ifndef GI_MOCK_WRITER_HPP
#define GI_MOCK_WRITER_HPP

#include "writers/writer.hpp"

namespace Spike
{
    struct MockWriter : public Writer
    {
        inline void close() {}
        inline void open(const std::string &file) {}
        inline void write(const std::string &line) override {}
    };
}

#endif