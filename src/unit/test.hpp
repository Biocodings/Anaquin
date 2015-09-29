#ifndef TEST_HPP
#define TEST_HPP

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

        static void fusionA();

        static void ladderA();
        static void ladderB();
        static void ladderAB();

        static void meta();
        
        static Test test(const std::string &);
    };
}

#endif