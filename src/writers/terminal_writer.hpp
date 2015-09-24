#ifndef TERMINAL_WRITER_HPP
#define TERMINAL_WRITER_HPP

#include <iostream>
#include "writers/writer.hpp"

namespace Anaquin
{
    class TerminalWriter : public Writer
    {
        public:
            inline void close() override
            {
                // Empty Implementation
            }
        
            inline void open(const FileName &file) override
            {
                // Empty Implementation
            }

            inline void write(const std::string &str) override
            {
                std::cout << str << std::endl;
            }
    };
}

#endif