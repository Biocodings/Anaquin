#include <vector>
#include <string>
#include <sstream>
#include "test.hpp"
#include "data/standard.hpp"
#include <boost/algorithm/string.hpp>

// Defined in main.cpp
extern int parse_options(int argc, char ** argv);

extern std::string RnaDataMixA();
extern std::string RnaDataMixB();
extern std::string RnaDataMixF();
extern std::string RnaDataMixAB();
extern std::string RnaStandGTF();

using namespace Anaquin;

void Test::clear()
{
    Standard::instance(true);
}

void Test::transA()
{
    Test::clear();
    Standard::instance().addRRef(Reader(RnaStandGTF(), DataMode::String));
    Standard::instance().addRMix(Reader(RnaDataMixA(), DataMode::String));
    Standard::instance().r_rna.finalize();
}

void Test::RnaQuin_B()
{
    Test::clear();
    Standard::instance().addRRef(Reader(RnaStandGTF(), DataMode::String));
    Standard::instance().addRMix(Reader(RnaDataMixB(), DataMode::String));
    Standard::instance().r_rna.finalize();
}

void Test::RnaQuin_AB()
{
    Test::clear();
    Standard::instance().addRRef(Reader(RnaStandGTF(), DataMode::String));
    Standard::instance().addRDMix(Reader(RnaDataMixAB(), DataMode::String));
    Standard::instance().r_rna.finalize();
}

void Test::RnaFoldChange()
{
    Test::clear();
    Standard::instance().addRDMix(Reader(RnaDataMixAB(), DataMode::String));
    Standard::instance().r_rna.finalize();
}

Test Test::test(const std::string &command)
{
    std::vector<std::string> tokens;
    boost::split(tokens, command, boost::is_any_of(" "));

    char *argv[tokens.size() + 1];

    argv[0] = new char(10);
    strcpy(argv[0], "anaquin");

    for (std::size_t i = 0; i < tokens.size(); i++)
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
    
    // Clear off anything made by a unit-test
    Standard::instance(true);

    return t;
}
