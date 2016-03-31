#ifndef SCRIPT_HPP
#define SCRIPT_HPP

#include <fstream>
#include <iostream>
#include <stdexcept>

// Defined in resources.cpp
extern std::string ViewerScript();

// Defined in resources.cpp
extern std::string ReportScript();

namespace Anaquin
{
    class Script
    {
        public:

            static void run(const std::string &code, const std::string &prefix, const std::string &args)
            {
                const auto f = [&](const std::string &cmd)
                {
                    const int status = system(cmd.c_str());
                    
                    if (status != 0)
                    {
                        throw std::runtime_error("Failed: " + cmd);
                    }
                };

                // Create a copy of the script
                const auto script = tmpnam(NULL);

                std::ofstream out(script);
                out << code;
                out.close();
                
                // Run the script with given arguments
                f(prefix + " " + script + " " + args);

                remove(script);
            }

            static void report(const std::string &args)
            {
                Script::run(ReportScript(), "python", args);            
            }
        
            static void viewer(const std::string &args)
            {
                Script::run(ViewerScript(), "python", args);
            }
    };
}

#endif