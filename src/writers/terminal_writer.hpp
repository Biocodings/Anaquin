#ifndef GI_TERMINAL_WRITER_HPP
#define GI_TERMINAL_WRITER_HPP

#include <iostream>
#include "writers/writer.hpp"

namespace Anaquin
{
    class TerminalWriter : public Writer
    {
        public:
            void close() override
            {
                // Empty Implementation
            }
        
            void open(const std::string &file) override
            {
                // Empty Implementation
            }

            void write(const std::string &str) override
            {
                std::cout << str << std::endl;
            }
    };
}

#endif