#ifndef GI_FILE_WRITER_HPP
#define GI_FILE_WRITER_HPP

#include <memory>
#include <fstream>
#include "writer.hpp"

namespace Spike
{
    class FileWriter : public Writer
    {
        public:
            FileWriter(const std::string &path, const std::string &file);

            void write(const std::string &line) override;

        private:
            std::shared_ptr<std::ofstream> _o;
    };
}

#endif