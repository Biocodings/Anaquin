#ifndef GI_TEST_HPP
#define GI_TEST_HPP

#include <string>
#include <catch.hpp>

namespace Anaquin
{
    struct Test
    {
        // Standard output
        std::string output;
        
        // Standard error
        std::string error;
        
        int status;

        static void transA();
        static void transB();
        static void transAB();

        static void variantA();
        static void variantF();

        static void meta();
        static void fusion();
        static void ladder();
        
        static Test test(const std::string &);
    };
}

#endif