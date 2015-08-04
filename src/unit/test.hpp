#ifndef GI_TEST_HPP
#define GI_TEST_HPP

#include <string>

namespace Anaquin
{
    struct CommandTest
    {
        // Standard output
        std::string output;
        
        // Standard error
        std::string error;
        
        int status;

        static CommandTest test(const std::string &);
    };
}

#endif