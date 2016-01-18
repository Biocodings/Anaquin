#ifndef TERMINAL_WRITER_HPP
#define TERMINAL_WRITER_HPP

#include <iostream>
#include "writers/writer.hpp"

namespace Anaquin
{
    class TerminalWriter : public Writer
    {
        public:

            inline void close() override {}

            inline void open(const FileName &) override {}

            inline void write(const std::string &str, bool) override
            {
                std::cout << str << std::endl;
            }

            inline void create(const std::string &) override {}
    };
}

#endif