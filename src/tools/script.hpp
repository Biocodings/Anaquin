#ifndef SCRIPT_TOOL_HPP
#define SCRIPT_TOOL_HPP

#include <string>
#include <sstream>
#include <ctype.h>

namespace Anaquin
{
    struct ScriptTool
    {
        static std::string clean(const std::string &str)
        {
            std::stringstream ss;
            
            for (auto i = 0; i < str.size(); i++)
            {
                auto curr = str.at(i);
                auto next = i != str.size()-1 ? str.at(i+1) : 0;

                if (curr == '$' && !isprint(next))
                {
                    continue;
                }
                
                ss << curr;
            }
            
            return ss.str();
        }        
    };
}

#endif