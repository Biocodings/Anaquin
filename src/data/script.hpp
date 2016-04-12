#ifndef SCRIPT_HPP
#define SCRIPT_HPP

#include <fstream>
#include <iostream>
#include <stdexcept>
#include "data/path.hpp"
#include <boost/format.hpp>

// Defined in resources.cpp
extern std::string ViewerScript();

// Defined in resources.cpp
extern std::string ReportScript();

// Defined in resources.cpp
extern std::string __anaquin__;

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
                const auto tmp = tmpFile();

                std::ofstream out(tmp);
                out << code;
                out.close();
                
                // Run the script with given arguments
                f(prefix + " " + tmp + " " + args);
            }

            template <typename Option> static void report(const std::string &type,
                                                          const std::string &file1,
                                                          const std::string &file2,
                                                          const Option &o)
            {
                const auto cmd = ((boost::format("%1% %2% %3% %4% %5% %6% %7%")
                                                    % type
                                                    % __anaquin__
                                                    % o.work
                                                    % o.mix
                                                    % o.index
                                                    % file1
                                                    % file2)).str();
                Script::run(ReportScript(), "python", cmd);
            }

            static void viewer(const std::string &args)
            {
                Script::run(ViewerScript(), "python", args);
            }
    };
}

#endif
