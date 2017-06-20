#include <mutex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "tools/system.hpp"

using namespace Anaquin;

bool System::isEmpty(const FileName& file)
{
    std::ifstream ss(file);
    return ss.peek() == std::ifstream::traits_type::eof();
}

FileName System::tmpFile()
{
    static std::mutex mtx;
    mtx.lock();
    
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wdeprecated-declarations"
    auto f = std::string(tmpnam(NULL));
    #pragma clang diagnostic pop

    auto random = [&](std::size_t len)
    {
        auto randchar = []() -> char
        {
            const char charset[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
            const size_t max_index = (sizeof(charset) - 1);
            return charset[rand() % max_index];
        };
        
        std::string str(len, 0);
        std::generate_n(str.begin(), len, randchar);
        
        return str;
    };

    // Avoid filename crashing...
    f = f + random(3) + random(3);
    
    mtx.unlock();
    
    return f;
}

void System::runCmd(const std::string &cmd)
{
    system(cmd.c_str());
}
