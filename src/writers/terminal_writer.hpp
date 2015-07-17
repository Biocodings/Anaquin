#ifndef GI_TERMINAL_WRITER_HPP
#define GI_TERMINAL_WRITER_HPP

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
        
            inline void open(const std::string &file) override
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