#ifndef GI_SCRIPT_HPP
#define GI_SCRIPT_HPP

#include <fstream>
#include <iostream>
#include <stdexcept>

extern std::string ViewerScript();

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
                const auto script = create(code);
                
                // Run the script with given arguments
                f(prefix + " " + script + " " + args);
            }

            static void viewer(const std::string &args)
            {
                return Script::run(ViewerScript(), "python", args);
            }

        private:

            static std::string create(const std::string &code)
            {
                std::ofstream out("python.py");
                out << code;
                out.close();

                return "python.py";
            }
    };
}

#endif