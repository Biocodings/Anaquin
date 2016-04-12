#ifndef PATH_HPP
#define PATH_HPP

#include <mutex>
#include <cstdio>
#include <algorithm>
#include "data/types.hpp"

namespace Anaquin
{
    inline std::string random(size_t len)
    {
        auto randchar = []() -> char
        {
            const char charset[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                   "abcdefghijklmnopqrstuvwxyz";
            const size_t max_index = (sizeof(charset) - 1);
            return charset[ rand() % max_index ];
        };

        std::string str(len, 0);
        std::generate_n(str.begin(), len, randchar);

        return str;
    }

    inline FileName tmpFile()
    {
        static std::mutex mtx;
        
        mtx.lock();
        
        auto f = std::string(tmpnam(NULL));

        // Avoid filename crashing...
        f = f + random(3) + random(3);
        
        mtx.unlock();
        
        return f;
    }
}

#endif