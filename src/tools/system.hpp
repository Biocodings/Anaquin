#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <string>
#include <sstream>
#include "data/data.hpp"
#include "tools/errors.hpp"

namespace Anaquin
{
    typedef std::string Package;
    typedef std::string Executable;
    
    inline FileName path2file(const Path &path)
    {
        auto tmp = path;
        const auto last = path.find_last_of("\\/");
        
        if (std::string::npos != last)
        {
            tmp.erase(0, last + 1);
        }
        
        return tmp;
    }
    
    inline std::vector<FileName> path2file(const std::vector<FileName> &files)
    {
        auto tmp = files;
        
        for (auto i = 0; i < tmp.size(); i++)
        {
            tmp[i] = path2file(tmp[i]);
        }

        return tmp;
    }
    
    struct System
    {
        // Check if the bash command is available
        static bool checkConsole(const Executable &);
        
        // Check if the R-package is available
        static bool checkRPack(const Package &);

        // Empty file?
        static bool isEmpty(const FileName& file);

        static void runCmd(const std::string &);
        static void runRScript(const std::string &code);
        
        static FileName tmpFile();

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
