#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "data/data.hpp"
#include "tools/errors.hpp"

namespace Anaquin
{
    typedef std::string Package;
    typedef std::string Executable;
    
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