#ifndef FILE_WRITER_HPP
#define FILE_WRITER_HPP

#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include "writers/writer.hpp"

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
                if (!path.empty())
                {
                    system((boost::format("mkdir -p %1%") % path).str().c_str());
                }
                
                _o = std::shared_ptr<std::ofstream>(new std::ofstream(!path.empty() ? path + "/" + file : file));
                
                if (!_o->good())
                {
                    throw std::runtime_error("Failed to open: " + file);
                }
            }

            inline void write(const std::string &line) override
            {
                *(_o) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << line << std::endl;
            }

            inline void create(const std::string &dir) override
            {
                const auto r = mkdir((path + "/" + dir).c_str(), 0666);

                if (r == -1 && errno != EEXIST)
                {
                    throw std::runtime_error("Failed to create directory: " + (path + "/" + dir) + " error: " + std::to_string(errno));
                }
            }
                
        private:
            std::string path;
            std::shared_ptr<std::ofstream> _o;
    };
}

#endif