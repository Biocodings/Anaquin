#include "writers/path_writer.hpp"

using namespace Spike;

void PathWriter::open(const std::string &file)
{
    _o = std::shared_ptr<std::ofstream>(new std::ofstream(path + "/" + file));
    
    if (!_o->good())
    {
        throw std::runtime_error("Failed to open: " + file);
    }
}

void PathWriter::write(const std::string &line)
{
    *(_o) << line << std::endl;
}