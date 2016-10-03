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

extern std::string VarDataBed();
extern std::string VarDataVCF();
extern std::string VarDataMixA();
extern std::string VarDataMixF();

using namespace Anaquin;

void Test::clear()
{
    Standard::instance(true);
}

void Test::transA()
{
    Test::clear();
    Standard::instance().addTRef(Reader(RnaStandGTF(), DataMode::String));
    Standard::instance().addTMix(Reader(RnaDataMixA(), DataMode::String));
    Standard::instance().r_rna.finalize();
}

void Test::RnaQuin_B()
{
    Test::clear();
    Standard::instance().addTRef(Reader(RnaStandGTF(), DataMode::String));
    Standard::instance().addTMix(Reader(RnaDataMixB(), DataMode::String));
    Standard::instance().r_rna.finalize();
}

void Test::RnaQuin_AB()
{
    Test::clear();
    Standard::instance().addTRef(Reader(RnaStandGTF(), DataMode::String));
    Standard::instance().addTDMix(Reader(RnaDataMixAB(), DataMode::String));
    Standard::instance().r_rna.finalize();
}

void Test::RnaFoldChange()
{
    Test::clear();
    Standard::instance().addTDMix(Reader(RnaDataMixAB(), DataMode::String));
    Standard::instance().r_rna.finalize();
}

void Test::VarQuinBed()
{
    Test::clear();
    Standard::instance().addVStd(Reader(VarDataBed(),  DataMode::String));
    Standard::instance().r_var.finalize();
}

void Test::variantA()
{
    Test::clear();
    Standard::instance().addVVar(Reader(VarDataVCF(),  DataMode::String));
    Standard::instance().addVStd(Reader(VarDataBed(),  DataMode::String));
    Standard::instance().addVMix(Reader(VarDataMixA(), DataMode::String));
    Standard::instance().r_var.finalize();
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
