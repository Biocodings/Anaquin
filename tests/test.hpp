#ifndef TEST_HPP
#define TEST_HPP

#include <string>
#include <catch.hpp>
#include "data/standard.hpp"
#include <boost/algorithm/string.hpp>

extern std::string RnaDataMixA();
extern std::string RnaDataMixAB();
extern std::string RnaStandGTF();

// Defined in main.cpp
extern int parse_options(int argc, char ** argv);

namespace Anaquin
{
    struct Test
    {
        // Standard output
        std::string output;
        
        // Standard error
        std::string error;
        
        int status;
    };

    inline void clrTest()
    {
        Standard::instance(true);
    }

    /*
     * RnaQuin analayis
     */
    
    inline void transA()
    {
        clrTest();
        Standard::instance().addRRef(Reader(RnaStandGTF(), DataMode::String));
        Standard::instance().addRMix(Reader(RnaDataMixA(), DataMode::String));
        Standard::instance().r_rna.finalize();
    }

    inline void RnaQuin_AB()
    {
        clrTest();
        Standard::instance().addRRef(Reader(RnaStandGTF(), DataMode::String));
        Standard::instance().addRDMix(Reader(RnaDataMixAB(), DataMode::String));
        Standard::instance().r_rna.finalize();
    }

    inline void RnaFoldChange()
    {
        clrTest();
        Standard::instance().addRDMix(Reader(RnaDataMixAB(), DataMode::String));
        Standard::instance().r_rna.finalize();        
    }
    
    inline Test runTest(const std::string &cmd)
    {
        clrTest();

        std::vector<std::string> tokens;
        boost::split(tokens, cmd, boost::is_any_of(" "));
        
        char *argv[tokens.size() + 1];
        
        argv[0] = new char(10);
        strcpy(argv[0], "anaquin");
        
        for (auto i = 0; i < tokens.size(); i++)
        {
            argv[i+1] = new char[tokens[i].size() + 1];
            strcpy(argv[i+1], tokens[i].c_str());
        }
        
        std::stringstream outputBuffer, errorBuffer;
        
        std::streambuf * const _errorBuffer  = std::cerr.rdbuf(errorBuffer.rdbuf());
        std::streambuf * const _outputBuffer = std::cout.rdbuf(outputBuffer.rdbuf());
        
        Test t;
        
        t.status = parse_options(static_cast<int>(tokens.size()) + 1, argv);
        t.error  = errorBuffer.str();
        t.output = outputBuffer.str();
        
        std::cout.rdbuf(_errorBuffer);
        std::cout.rdbuf(_outputBuffer);
        
        for (std::size_t i = 0; i < tokens.size(); i++)
        {
            delete[] argv[i];
        }
        
        clrTest();
        
        return t;
    }
}

#endif
