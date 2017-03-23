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

        static void clear();

        /*
         * RnaQuin analayis
         */
        
        static void transA();
        static void RnaQuin_B();
        static void RnaQuin_AB();
        static void RnaFoldChange();

        /*
         * Variant analayis
         */

        /*
         * Fusion analayis
         */
        
        /*
         * Ladder analayis
         */
        
        /*
         * Metagenomics analayis
         */
        
        static Test test(const std::string &);
    };
}

#endif
