#ifndef GI_PATH_WRITER_HPP
#define GI_PATH_WRITER_HPP

#include <memory>
#include <fstream>
#include "writer.hpp"

namespace Spike
{
    class PathWriter : public Writer
    {
        public:
            PathWriter(const std::string &path) : path(path) {}

            inline void close()
            {
                _o->close();
                _o.reset();
            }

            void open(const std::string &file) override;
            void write(const std::string &line) override;

        private:
            std::string path;
            std::shared_ptr<std::ofstream> _o;
    };
}

#endif