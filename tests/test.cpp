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

extern std::string FusionDataMixA();
extern std::string FusionDataRef();

extern std::string TransDataMixA();
extern std::string TransDataMixB();
extern std::string TransDataMixAB();
extern std::string TransStandGTF();

extern std::string MetaDataBed();
extern std::string MetaDataMix();

extern std::string VarDataBed();
extern std::string VarDataVCF();
extern std::string VarDataMixA();
extern std::string VarDataMixF();

using namespace Anaquin;

void Test::fusionA()
{
    Standard::instance(true);
    Standard::instance().f_ref(Reader(FusionDataRef(), DataMode::String));
    Standard::instance().f_mix(Reader(FusionDataMixA(), DataMode::String));
    Standard::instance().r_fus.validate();
}

void Test::meta()
{
    Standard::instance(true);
    Standard::instance().m_mix_1(Reader(MetaDataMix(), DataMode::String));
    Standard::instance().m_mix_2(Reader(MetaDataMix(), DataMode::String));
    Standard::instance().r_meta.validate();
}

void Test::ladder()
{
    Standard::instance(true);
    Standard::instance().l_mix(Reader(LadderDataMix(), DataMode::String));
    Standard::instance().r_lad.validate();
}

void Test::transA()
{
    Standard::instance(true);
    Standard::instance().r_ref(Reader(TransStandGTF(), DataMode::String));
    Standard::instance().r_mix(Reader(TransDataMixA(), DataMode::String));
    Standard::instance().r_trans.validate();
}

void Test::transB()
{
    Standard::instance(true);
    Standard::instance().r_ref(Reader(TransStandGTF(), DataMode::String));
    Standard::instance().r_mix(Reader(TransDataMixB(), DataMode::String));
    Standard::instance().r_trans.validate();
}

void Test::transAB()
{
    Standard::instance(true);
    Standard::instance().r_ref(Reader(TransStandGTF(),  DataMode::String));
    Standard::instance().r_mix(Reader(TransDataMixAB(), DataMode::String));
    Standard::instance().r_trans.validate();
}

void Test::variantA()
{
    Standard::instance(true);
    Standard::instance().v_var(Reader(VarDataBed(),  DataMode::String));
    //Standard::instance().v_std(Reader(VarStandGTF(), DataMode::String));
    Standard::instance().v_mix(Reader(VarDataMixA(),  DataMode::String));
    Standard::instance().r_var.validate();
}

void Test::variantF()
{
    Standard::instance(true);
    Standard::instance().v_var(Reader(VarDataBed(),  DataMode::String));
    //Standard::instance().v_std(Reader(VarStandGTF(), DataMode::String));
    Standard::instance().v_mix(Reader(VarDataMixA(),  DataMode::String));
    Standard::instance().r_var.validate();
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