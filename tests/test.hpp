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

        static void VarQuinBed();
        static void variantA();
        static void variantF();

        /*
         * Fusion analayis
         */
        
        static void fusionA();

        /*
         * Ladder analayis
         */
        
        static void ladderA();
        static void ladderB();
        static void ladderAB();

        /*
         * Metagenomics analayis
         */
        
        static void meta();
        
        static Test test(const std::string &);
    };
}

#endif
