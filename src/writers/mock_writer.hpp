#ifndef GI_MOCK_WRITER_HPP
#define GI_MOCK_WRITER_HPP

#include <vector>
#include "writers/writer.hpp"

namespace Spike
{
    struct MockWriter : public Writer
    {
        inline void write(const std::string &line) override
        {
            lines.push_back(line);
        }

        std::vector<std::string> lines;
    };
}

#endif