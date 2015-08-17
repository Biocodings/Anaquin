#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "unit/test.hpp"
#include "data/standard.hpp"
#include <boost/algorithm/string.hpp>

// Defined in main.cpp
extern int parse_options(int argc, char ** argv);

extern std::string LadderDataMix();

extern std::string FusionDataMix();
extern std::string FusionDataRef();
extern std::string FusionNormalRef();
extern std::string FusionMutatedRef();

extern std::string TransDataMix();
extern std::string TransStandGTF();

extern std::string MetaDataFA();
extern std::string MetaDataBed();
extern std::string MetaDataMix();
extern std::string MetaDataTab();

extern std::string VarDataBed();
extern std::string VarDataMix();
extern std::string VarStandGTF();

using namespace Anaquin;

void Test::fusion()
{
    Standard::instance().f_ref(Reader(FusionDataRef(), DataMode::String));
    Standard::instance().f_mix(Reader(FusionDataMix(), DataMode::String));
}

void Test::meta()
{
    Standard::instance().m_ref(Reader(MetaDataBed(), DataMode::String));
    Standard::instance().m_mix_1(Reader(MetaDataMix(), DataMode::String));
    Standard::instance().m_mix_2(Reader(MetaDataMix(), DataMode::String));
}

void Test::ladder()
{
    Standard::instance().l_mix(Reader(LadderDataMix(), DataMode::String));
}

void Test::trans()
{
    Standard::instance().r_ref(Reader(TransStandGTF(), DataMode::String));
    Standard::instance().r_mix(Reader(TransDataMix(), DataMode::String));
}

void Test::variant()
{
    Standard::instance().v_var(Reader(VarDataBed(),  DataMode::String));
    Standard::instance().v_std(Reader(VarStandGTF(), DataMode::String));
    Standard::instance().v_mix(Reader(VarDataMix(),  DataMode::String));
}

Test Test::test(const std::string &command)
{
    std::vector<std::string> tokens;
    boost::split(tokens, command, boost::is_any_of(" "));

    char *argv[tokens.size() + 1];

    argv[0] = new char(10);
    strcpy(argv[0], "test.exe");

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