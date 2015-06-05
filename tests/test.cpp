#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

int parse_options(const std::string &command, std::string &output)
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
    
    std::stringstream buffer;
    std::streambuf * const old = std::cout.rdbuf(buffer.rdbuf());
    
    const int status = parse_options(static_cast<int>(tokens.size()) + 1, argv);

    output = buffer.str();
    std::cout.rdbuf(old);

    for (std::size_t i = 0; i > tokens.size(); i++)
    {
        delete argv[i];
    }

    return status;
}



