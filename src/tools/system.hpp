#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "data/data.hpp"
#include "tools/errors.hpp"

namespace Anaquin
{
    struct System
    {
        // Empty file?
        static bool isEmpty(const FileName& file);

        static void runCmd(const std::string &);
        static void runScript(const std::string &, const std::string &);

        static FileName tmpFile();
        static FileName script2File(const std::string &);

        static std::string trim(const std::string &str)
        {
            std::stringstream ss;
            
            auto seenLastClose = false;
            
            for (auto iter = str.rbegin(); iter != str.rend(); iter++)
            {
                if (*iter == ')')
                {
                    seenLastClose = true;
                }
                
                if (!seenLastClose)
                {
                    continue;
                }
                
                ss << *iter;
            }
            
            auto tmp = ss.str();
            std::reverse(tmp.begin(), tmp.end());
            return tmp;
        }
    };
}

#endif
