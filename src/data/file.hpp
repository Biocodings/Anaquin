#ifndef FILE_HPP
#define FILE_HPP

#include <mutex>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "data/data.hpp"

namespace Anaquin
{
    inline bool isEmpty(const FileName& file)
    {
        std::ifstream ss(file);
        return ss.peek() == std::ifstream::traits_type::eof();
    }

    inline std::string random(size_t len)
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
    }

    inline FileName tmpFile()
    {
        static std::mutex mtx;        
        mtx.lock();
        
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
        auto f = std::string(tmpnam(NULL));
#pragma clang diagnostic pop

        // Avoid filename crashing...
        f = f + random(3) + random(3);
        
        mtx.unlock();
        
        return f;
    }
}

#endif