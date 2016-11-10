#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
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
        static void runScript(const std::string &, const std::string &, const std::string &) throw(FailedCommandException);

        static FileName tmpFile();
    };
}

#endif
