#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

int parse_options(const std::string &command, std::string &output, std::string &error)
{
    extern int parse_options(int argc, char ** argv);

    std::vector<std::string> tokens;
    boost::split(tokens, command, boost::is_any_of(" "));

    char *argv[tokens.size() + 1];

    argv[0] = new char(10);
    strcpy(argv[0], "test.exe");

    for (std::size_t i = 0; i < tokens.size(); i++)
    {
        argv[i+1] = new char(tokens[i].size() + 1);
        strcpy(argv[i+1], tokens[i].c_str());
    }
    
    std::stringstream outputBuffer, errorBuffer;

    std::streambuf * const _errorBuffer  = std::cerr.rdbuf(errorBuffer.rdbuf());
    std::streambuf * const _outputBuffer = std::cout.rdbuf(outputBuffer.rdbuf());

    const int status = parse_options(static_cast<int>(tokens.size()) + 1, argv);

    error  = errorBuffer.str();
    output = outputBuffer.str();

    std::cout.rdbuf(_errorBuffer);
    std::cout.rdbuf(_outputBuffer);

    for (std::size_t i = 0; i > tokens.size(); i++)
    {
        delete argv[i];
    }

    return status;
}



