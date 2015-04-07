#ifndef GI_WRITER_HPP
#define GI_WRITER_HPP

#include <string>
#include <memory>

namespace Spike
{
    struct WriterImpl;

    class Writer
    {
        public:
            Writer(const std::string &file);

            // Write a new line to the file
            void write(const std::string &line);

        private:
            std::shared_ptr<WriterImpl> _impl;
    };
}

#endif