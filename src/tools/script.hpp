#ifndef SCRIPT_HPP
#define SCRIPT_HPP

#include <string>
#include <sstream>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "tools/errors.hpp"
#include "tools/system.hpp"
#include <boost/format.hpp>

// Defined in resources.cpp
extern std::string ReportScript();

// Defined in resources.cpp
extern std::string __anaquin__;

namespace Anaquin
{
    struct Script
    {
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