#ifndef FILE_WRITER_HPP
#define FILE_WRITER_HPP

#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include "tools/script.hpp"
#include "writers/writer.hpp"
#include <boost/algorithm/string/predicate.hpp>

namespace Anaquin
{
    class FileWriter : public Writer
    {
        public:

            FileWriter(const std::string &path) : path(path) {}

            inline void close() override
            {
                _o->close();
                _o.reset();
            }

            inline void open(const FileName &file) override
            {
                isScript = boost::algorithm::ends_with(file, ".R") ||
                           boost::algorithm::ends_with(file, ".py");

                if (!path.empty())
                {
                    system((boost::format("mkdir -p %1%") % path).str().c_str());
                }
                
                const auto target = !path.empty() ? path + "/" + file : file;
                _o = std::shared_ptr<std::ofstream>(new std::ofstream(target));
                
                if (!_o->good())
                {
                    throw std::runtime_error("Failed to open: " + target);
                }
            }

            inline void write(const std::string &x, bool newLine = true) override
            {
                if (isScript)
                {
                    *(_o) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << Script::trim(x);
                }
                else
                {
                    *(_o) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << x;
                }

                if (newLine)
                {
                    *(_o) << std::endl;
                }
            }

            inline void create(const std::string &dir) override
            {
                const auto r = mkdir((path + "/" + dir).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

                if (r == -1 && errno != EEXIST)
                {
                    throw std::runtime_error("Failed to create directory: " + (path + "/" + dir) + " error: " + std::to_string(errno));
                }
            }
                
        private:
        
            // Trim script?
            bool isScript;
        
            std::string path;
            std::shared_ptr<std::ofstream> _o;
    };
}

#endif