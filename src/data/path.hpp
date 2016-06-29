#ifndef PATH_HPP
#define PATH_HPP

#include <mutex>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "data/data.hpp"

namespace Anaquin
{
    inline FileName path2file(const FileName &file)
    {
        auto tmp = file;
        
        // Remove directory if present.
        const size_t last_slash_idx = tmp.find_last_of("\\/");
        
        if (std::string::npos != last_slash_idx)
        {
            tmp.erase(0, last_slash_idx + 1);
        }

        return tmp;
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