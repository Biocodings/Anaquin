#ifndef GI_TEST_HPP
#define GI_TEST_HPP

#include <string>

namespace Anaquin
{
    struct Test
    {
        // Standard output
        std::string output;
        
        // Standard error
        std::string error;
        
        int status;

        // Apply default resources to fusion
        static void fusion();

        static Test test(const std::string &);
    };
}

#endif