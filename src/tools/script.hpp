#ifndef SCRIPT_HPP
#define SCRIPT_HPP

#include <string>
#include <sstream>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "data/file.hpp"
#include <boost/format.hpp>

// Defined in resources.cpp
extern std::string ReportScript();

// Defined in resources.cpp
extern std::string __anaquin__;

namespace Anaquin
{
    struct Script
    {
        static void run(const std::string &code, const std::string &prefix, const std::string &args, const std::string &exts = "")
        {
            const auto f = [&](const std::string &cmd)
            {
                std::cout << cmd << std::endl;
                
                const int status = system(cmd.c_str());
                
                if (status != 0)
                {
                    throw std::runtime_error("Failed: " + cmd);
                }
            };
            
            // Create a copy of the script
            const auto tmp = tmpFile() + exts;
            
            std::ofstream out(tmp);
            out << code;
            out.close();
            
            // Run the script with given arguments
            f(prefix + " " + tmp + " " + args);
        }
        
        // Eg: python reports.py RnaQuin reports.pdf output
        template <typename Option> static void report(const std::string &type,
                                                      const FileName &report,
                                                      const FileName &inputs,
                                                      const Option &o)
        {
            const auto cmd = ((boost::format("%1% %2% %3%") % type
                                                            % (o.work + "/" + report)
                                                            % inputs)).str();
            Script::run(ReportScript(), "python", cmd);
        }
        
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